import glob
import json
from collections import defaultdict

from .utils import *


class GenomeProfiler:
    '''
    Profile taxonomic genomes using a set of marker genes.
    '''
    def __init__(self, file, output, threads=32):
        self.file = file
        self.output = output
        self.threads = threads

        self.aset = {'l2', 'l11', 'l10e', 'l15e', 'l18e', 's3ae', 's19e', 's28e'}
        self.bset = {'l2', 'l11', 'l20', 'l27', 's2', 's7', 's9', 's16'}
        self.nset = set()

        self.ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

        ## genome copies
        self.copies = {'bacteria': 0, 'archaea': 0}

        ## temporary variables for diamond's hits or minimap2's mappings
        self.hits, self.maps = [], []

        ## taxonomy assignments
        self.assignments = {}


    def run_kraken(self, db_kraken):
        '''
        Run kraken2 for pre-filtering.
        '''
        subprocess.run([
            'kraken2',
            '--db', db_kraken,
            '--report', get_filename(self.file, self.output, '.kraken.report.tmp'),
            '--output', get_filename(self.file, self.output, '.kraken.output.tmp'),
            '--threads', str(self.threads),
            self.file,
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    def parse_kraken(self):
        '''
        Parse kraken2's output by
            1: recording the ids of eukaryota (2759), viruses (10239), other entries (2787854).
            2: adding the id of the sequence to a negative set.
        '''
        seqid, taxid = set(), set()
        record = False
        with open(get_filename(self.file, self.output, '.kraken.report.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                if ls[3] in {'D', 'R1'}:
                    record = ls[5].strip() in {'Eukaryota', 'Viruses', 'other entries'}

                if record:
                    taxid.add(ls[4])

        with open(get_filename(self.file, self.output, '.kraken.output.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                if ls[2] in taxid:
                    seqid.add(ls[1])

        self.nset.update(seqid)


    def run_diamond(self, db, max_target_seqs=25, evalue=1e-15, identity=0, subject_cover=75):
        '''
        Run diamond to get total prokaryotic genome copies.
        '''
        outfmt = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send', 'evalue', 'bitscore']
        subprocess.run([
            'diamond', 'blastx',
            '--db', os.path.join(db, 'prot.dmnd'),
            '--query', self.file,
            '--out', get_filename(self.file, self.output, '.diamond.tmp'),
            '--outfmt', '6', *outfmt,
            '--evalue', str(evalue), '--subject-cover', str(subject_cover), '--id', str(identity),
            '--range-culling', '-F', '15', '--range-cover', '25',
            '--max-hsps', '0', '--max-target-seqs', str(max_target_seqs),
            '--threads', str(self.threads)
        ], check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


    def parse_diamond(self):
        '''
        Parse diamond's output and record the hits.
        '''
        qrange = defaultdict(set)
        srange = {}

        with open(get_filename(self.file, self.output, '.diamond.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qseqid, sseqid = ls[0], ls[1]
                qstart, qend = sort_coordinate(int(ls[5]), int(ls[6]))
                sstart, send = sort_coordinate(int(ls[8]), int(ls[9]))
                slen = int(ls[7])

                ## bypass non-prokaryotic reads
                if qseqid not in self.nset:
                    if (
                        qseqid not in qrange or
                        all(compute_overlap((qstart, qend, *x), max) < 0.25 for x in qrange[qseqid])
                    ):
                        qrange[qseqid].add((qstart, qend))

                        ss = sseqid.split('-')
                        gene, kingdom = ss[0], ss[-2]
                        if (
                            (gene in self.bset and kingdom == 'bacteria') or
                            (gene in self.aset and kingdom == 'archaea')
                        ):
                            if sseqid not in srange:
                                srange[sseqid] = np.zeros(slen)
                            srange[sseqid][range(sstart, send)] += 1

                            ## append qseqid and coordinates for back-tracing
                            self.hits.append([qseqid, kingdom, gene, qstart, qend])

        for key, val in srange.items():
            cut = round(len(val) / 4) # in this case same as counting hits since all hits have subject cover > 75%
            kingdom = key.split('-')[-2]
            if kingdom == 'bacteria':
                self.copies[kingdom] += np.mean(val[cut:-cut]) / len(self.bset)
            else:
                self.copies[kingdom] += np.mean(val[cut:-cut]) / len(self.aset)


    def run_minimap(self, db, secondary_num=2147483647, secondary_ratio=0.9):
        '''
        Run minimap2 to get taxonomic profiles.
        '''
        with open(get_filename(self.file, self.output, '.sequence.tmp'), 'w') as w:
            w.write(extract_sequences(self.file, {x[0] for x in self.hits}))

        ## consider each kingdom + gene separately
        genes = defaultdict(set)
        for i in self.hits:
            genes[i[1] + '.' + i[2].replace('/', '_')].add(i[0])

        with open(get_filename(self.file, self.output, '.minimap.tmp'), 'w') as w:
            for key, val in genes.items():
                sequences = extract_sequences(get_filename(self.file, self.output, '.sequence.tmp'), val)
                subprocess.run([
                    'minimap2',
                    '-cx', 'map-ont',
                    '-f', '0',
                    '-N', str(secondary_num), '-p', str(secondary_ratio),
                    '-t', str(self.threads),
                    os.path.join(db, 'nucl.' + key + '.mmi'), '-',
                ], check=True, stdout=w, stderr=subprocess.DEVNULL, input=sequences, text=True)


    def parse_minimap(self):
        '''
        Parse minimap2's output and record the mappings.
        '''
        coordinates = defaultdict(set)
        for i in self.hits:
            coordinates[i[0]].add(tuple(i[-2:]))

        with open(get_filename(self.file, self.output, '.minimap.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qstart, qend, qseqid, sseqid = int(ls[2]), int(ls[3]), ls[0], ls[5]

                AS = int(ls[14].split('AS:i:')[-1])
                MS = int(ls[13].split('ms:i:')[-1])
                ID = 1 - float(ls[19].split('de:f:')[-1]) if ls[16] == 'tp:A:S' or ls[16] == 'tp:A:i' else 1 - float(ls[20].split('de:f:')[-1])

                ## filter out non-overlapping mappings
                for i in coordinates[qseqid]:
                    if compute_overlap((qstart, qend, *i)) > 0:
                        self.maps.append([qseqid, sseqid, AS, MS, ID])


    def postprocess(self, db):
        '''
        Post-processing and label reassignment using EM.
        '''
        accession2lineage = {}
        with open(os.path.join(db, 'metadata.tsv')) as f:
            next(f)
            for line in f:
                ls = line.rstrip().split('\t')
                accession2lineage[ls[0]] = ';'.join(ls[1:])

        ## paste assigned taxonomy then sort
        maps = sorted([
            [*x, accession2lineage[x[1].rsplit('_', 1)[0]]] for x in self.maps
        ], key=lambda x: (x[0], x[2], x[3], x[4]), reverse=True)

        ## keep only the first per qseqid and lineage, remove all inferior alignments
        data = []
        duplicates = set()
        max_scores = defaultdict(lambda: {'AS': 0, 'MS': 0, 'ID': 0})

        for row in maps:
            max_scores[row[0]]['AS'] = max(max_scores[row[0]]['AS'], row[2])
            max_scores[row[0]]['MS'] = max(max_scores[row[0]]['MS'], row[3])
            max_scores[row[0]]['ID'] = max(max_scores[row[0]]['ID'], row[4])

            if (row[0], row[-1]) not in duplicates:
                data.append(row)
                duplicates.add((row[0], row[-1]))

        data = [row for row in data if (
            row[2] > max_scores[row[0]]['AS'] * 0.99 or
            row[3] > max_scores[row[0]]['MS'] * 0.99 or
            row[4] > max_scores[row[0]]['ID'] * 0.999
        )]

        ## create a matrix then fill
        qseqids, lineages = np.unique([row[0] for row in data]), np.unique([row[-1] for row in data])
        qseqid2index, lineage2index = {qseqid: index for index, qseqid in enumerate(qseqids)}, {lineage: index for index, lineage in enumerate(lineages)}

        matrix = np.zeros((len(qseqids), len(lineages)), dtype=int)
        for row in data:
            matrix[qseqid2index[row[0]], lineage2index[row[-1]]] += 1

        ## run EM using the count matrix as input
        assignments = reassign_taxonomy(matrix)
        ties = defaultdict(list)
        for qseqid, lineage in enumerate(assignments):
            if len(assignment := lineages[lineage]) > 1:
                ties[tuple(assignment)].append(qseqids[qseqid])
            else:
                self.assignments[qseqids[qseqid]] = assignment[0]

        ## resolve ties for equal probability cases using AS and ID
        for key, val in ties.items():
            target = [row for row in data if row[-1] in key and row[0] in val]

            scores = defaultdict(lambda: defaultdict(list))
            for row in target:
                scores[row[-1]]['AS'].append(row[2])
                scores[row[-1]]['MS'].append(row[3])
                scores[row[-1]]['ID'].append(row[4])

            ## if still tie in AS and de, choose the one with known species name
            target = sorted([
                [np.mean(val['AS']), np.mean(val['MS']), np.mean(val['ID']), not bool(re.search(' sp\.$| sp\. |sp[0-9]+', key.split(';')[-1])), key] for key, val in scores.items()
            ], reverse=True)[0][-1]

            for qseqid in val:
                self.assignments[qseqid] = target


    def run(self, db, db_kraken=None, skip_profile=False, skip_clean=False,
            max_target_seqs=25, evalue=1e-15, identity=0, subject_cover=75,
            secondary_num=2147483647, secondary_ratio=0.9):
        '''
        Run the pipeline.
        '''
        if db_kraken is not None:
            logger.info('Filtering reads ...')
            self.run_kraken(db_kraken)
            self.parse_kraken()
            logger.info('... removed {} putatively non-prokaryotic reads.'.format(len(self.nset)))

        logger.info('Estimating genome copies ...')
        self.run_diamond(db, max_target_seqs, evalue, identity, subject_cover)
        self.parse_diamond()
        logger.info('... found {} copies of genomes (bacteria: {}; archaea: {}).'.format(
            sum(self.copies.values()), self.copies['bacteria'], self.copies['archaea']))

        if not skip_profile:
            logger.info('Assigning taxonomy ...')
            self.run_minimap(db, secondary_num, secondary_ratio)
            self.parse_minimap()
            self.postprocess(db)

            ## fill missing ones according to hits
            replacement = {
                'bacteria': ';'.join(['2|Bacteria'] + ['0|unclassified Bacteria ' + x.lower() for x in self.ranks[1:]]),
                'archaea': ';'.join(['2157|Archaea'] + ['0|unclassified Archaea ' + x.lower() for x in self.ranks[1:]])
            }
            
            ## fit gtdb style
            if self.assignments and '|' not in next(iter(self.assignments.values())).split(';')[0]:
                replacement = {key: ';'.join(x.split('|')[-1] for x in val.split(';')) for key, val in replacement.items()}

            ## count assigned taxonomic labels
            counts, total_counts = defaultdict(lambda: 0), defaultdict(lambda: 0)
            for i in self.hits:
                counts[(self.assignments.get(i[0], replacement.get(i[1])), i[1])] += 1
                total_counts[i[1]] += 1

            # generate a profile output
            self.profile = sorted([
                [*key[0].split(';'), val * self.copies[key[1]] / total_counts[key[1]], val / sum(total_counts.values())] for key, val in counts.items()
            ], key=lambda x: (x[-2], x[-3]))

            richness = {'bacteria': 0, 'archaea': 0}
            with open(get_filename(self.file, self.output, '.tsv'), 'w') as w:
                w.write('\t'.join(self.ranks + ['copies', 'abundance']) + '\n')
                for line in self.profile:
                    if not re.search('unclassified (Bacteria|Archaea) species', line[6]):
                        richness[line[0].split('|')[-1].lower()] += 1
                    w.write('\t'.join(str(x) for x in line) + '\n')

            logger.info('... found {} unique species (bacteria: {}; archaea: {}).'.format(
                sum(richness.values()), richness['bacteria'], richness['archaea']))

        else:
            self.profile = sorted([
                ['2|Bacteria', self.copies['bacteria'], self.copies['bacteria'] / sum(self.copies.values())],
                ['2157|Archaea', self.copies['archaea'], self.copies['archaea'] / sum(self.copies.values())]
            ], key=lambda x: (x[-2], x[-3]))

            with open(get_filename(self.file, self.output, '.tsv'), 'w') as w:
                w.write('\t'.join(['superkingdom', 'copies', 'abundance']) + '\n')
                for line in self.profile:
                    if line[1] != 0:
                        w.write('\t'.join(str(x) for x in line) + '\n')

        ## save reads
        reads = {x:{'remark': 'putatively non-prokaryotic', 'lineage': 'others'} for x in self.nset}
        if not skip_profile:
            reads.update({x[0]:{'remark': 'marker-gene-containing', 'lineage': self.assignments.get(x[0], replacement.get(x[1]))} for x in self.hits})
        else:
            reads.update({x[0]:{'remark': 'marker-gene-containing', 'lineage': x[1]} for x in self.hits})

        with open(get_filename(self.file, self.output, '.json'), 'w') as w:
            json.dump(dict(sorted(reads.items())), w, indent=4)

        ## clean up
        if not skip_clean:
            for f in glob.glob(get_filename(self.file, self.output, '.*.tmp')):
                os.remove(f)
