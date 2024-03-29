import glob
import json
import numpy as np

from collections import defaultdict
from scipy.sparse import csr_matrix
from .utils import *


class GenomeProfiler:
    '''
    Profile taxonomic genomes using a set of marker genes.
    '''
    def __init__(self, file, db, output, threads=os.cpu_count()):
        self.file = file
        self.db = db
        self.output = output
        self.threads = threads

        self.aset = {'l2', 'l11', 'l10e', 'l15e', 'l18e', 's3ae', 's19e', 's28e'}
        self.bset = {'l2', 'l11', 'l20', 'l27', 's2', 's7', 's9', 's16'}
        self.nset = set()

        self.ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

        ## genome copies
        self.copies = {'bacteria': 0, 'archaea': 0}

        ## temporary variables for diamond's hits or minimap2's alignments
        self.hits, self.alignments = [], []

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
        qseqids, taxids = set(), set()
        record = False
        with open(get_filename(self.file, self.output, '.kraken.report.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                if ls[3] in {'D', 'R1'}:
                    record = ls[5].strip() in {'Eukaryota', 'Viruses', 'other entries'}

                if record:
                    taxids.add(ls[4])

        with open(get_filename(self.file, self.output, '.kraken.output.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                if ls[2] in taxids:
                    qseqids.add(ls[1])

        self.nset.update(qseqids)

    def run_diamond(self, max_target_seqs=25, evalue=1e-15, identity=0, subject_cover=75):
        '''
        Run diamond to get total prokaryotic genome copies.
        '''
        outfmt = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send', 'evalue']
        subprocess.run([
            'diamond', 'blastx',
            '--db', os.path.join(self.db, 'prot.dmnd'),
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
        qcoords = defaultdict(set)

        with open(get_filename(self.file, self.output, '.diamond.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qseqid, sseqid = ls[0], ls[1]
                qstart, qend = sort_coordinate(int(ls[5]), int(ls[6]))

                ## bypass non-prokaryotic reads
                if qseqid not in self.nset:
                    if (
                        qseqid not in qcoords or
                        all(compute_overlap((qstart, qend, *qcoord), max) < 0.25 for qcoord in qcoords[qseqid])
                    ):
                        qcoords[qseqid].add((qstart, qend))

                        ss = sseqid.split('-')
                        family, kingdom = ss[0], ss[-2]
                        if (
                            (family in self.bset and kingdom == 'bacteria') or
                            (family in self.aset and kingdom == 'archaea')
                        ):
                            ## append qseqid and coordinates for back-tracing
                            self.hits.append([qseqid, kingdom, family, qstart, qend])
                            self.copies[sseqid.split('-')[-2]] += 0.125

    def run_minimap(self, secondary_num=2147483647, secondary_ratio=0.9):
        '''
        Run minimap2 to get taxonomic profiles.
        '''
        with open(get_filename(self.file, self.output, '.sequence.tmp'), 'w') as w:
            w.write(extract_sequences(self.file, {hit[0] for hit in self.hits}))

        ## consider each kingdom + family separately
        qseqids = defaultdict(set)
        for hit in self.hits:
            qseqids[hit[1] + '.' + hit[2].replace('/', '_')].add(hit[0])

        with open(get_filename(self.file, self.output, '.minimap.tmp'), 'w') as w:
            for family, qseqid in qseqids.items():
                sequences = extract_sequences(get_filename(self.file, self.output, '.sequence.tmp'), qseqid)
                subprocess.run([
                    'minimap2',
                    '-cx', 'map-ont',
                    '-f', '0',
                    '-N', str(secondary_num), '-p', str(secondary_ratio),
                    '-t', str(self.threads),
                    os.path.join(self.db, 'nucl.' + family + '.mmi'), '-',
                ], check=True, stdout=w, stderr=subprocess.DEVNULL, input=sequences, text=True)

    def parse_minimap(self):
        '''
        Parse minimap2's output and record alignments.
        '''
        qcoords = defaultdict(set)
        for hit in self.hits:
            qcoords[hit[0]].add(tuple(hit[-2:]))

        with open(get_filename(self.file, self.output, '.minimap.tmp')) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qstart, qend, qseqid, sseqid = int(ls[2]), int(ls[3]), ls[0], ls[5]

                AS = int(ls[14].split('AS:i:')[-1])
                DE = 1 - float((ls[19] if ls[16] in {'tp:A:S', 'tp:A:i'} else ls[20]).split('de:f:')[-1]) # gap-compressed identity
                ID = int(ls[9]) / int(ls[10]) # gap-uncompressed identity

                ## filter out non-overlapping alignments
                for qcoord in qcoords[qseqid]:
                    if compute_overlap((qstart, qend, *qcoord)) > 0:
                        self.alignments.append([qseqid, sseqid, AS, DE, ID])

    def postprocess(self, max_iteration=1000, epsilon=1e-10):
        '''
        Post-processing and label reassignment using EM.
        '''
        accession2lineage = {}
        with open(os.path.join(self.db, 'metadata.tsv')) as f:
            next(f)
            for line in f:
                ls = line.rstrip().split('\t')
                accession2lineage[ls[0]] = ';'.join(ls[1:])

        ## alignment filtering based on AS, DE, and ID, keep only the first per qseqid and lineage, remove all inferior ones
        data = []
        duplicates = set()
        max_scores = defaultdict(lambda: {'AS': 0, 'DE': 0, 'ID': 0})

        for alignment in self.alignments:
            max_scores[alignment[0]]['AS'] = max(max_scores[alignment[0]]['AS'], alignment[2])
            max_scores[alignment[0]]['DE'] = max(max_scores[alignment[0]]['DE'], alignment[3])
            max_scores[alignment[0]]['ID'] = max(max_scores[alignment[0]]['ID'], alignment[4])

        for row in sorted(self.alignments, key=lambda alignment: (alignment[0], alignment[2], alignment[3], alignment[4]), reverse=True):
            if (
                max(row[2] / 0.9975, row[2] + 25) > max_scores[row[0]]['AS'] or
                row[3] / 0.9995 > max_scores[row[0]]['DE'] or
                row[4] / 0.9995 > max_scores[row[0]]['ID']
            ):
                species = accession2lineage[row[1].rsplit('_', 1)[0]]
                if (row[0], species) not in duplicates:
                    data.append(row + [species])
                    duplicates.add((row[0], species))

        # ## save pairwise gap-compressed identity for ANI calculation
        self.identities = {(row[0], row[-1]): (row[3], row[4]) for row in data}

        ## create a matrix then fill
        qseqids, lineages = np.unique([row[0] for row in data]), np.unique([row[-1] for row in data])
        qseqid2index = {qseqid: index for index, qseqid in enumerate(qseqids)}
        lineage2index = {lineage: index for index, lineage in enumerate(lineages)}

        rows = [qseqid2index[row[0]] for row in data]
        cols = [lineage2index[row[-1]] for row in data]
        matrix = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(qseqids), len(lineages)), dtype=int)

        ## run EM using the count matrix as input
        n_reads, n_mappings = matrix.shape

        ## init
        p_mappings = np.ones((1, n_mappings)) / n_mappings
        p_mappings_hist = p_mappings.copy()

        iteration = 0
        while iteration < max_iteration:
            iteration += 1

            ## e-step
            p_reads = matrix.multiply(p_mappings) / matrix.dot(p_mappings.T)

            ## m-step
            p_mappings = np.sum(p_reads, axis=0) / n_reads

            ## check convergence
            if np.sum(np.abs(p_mappings - p_mappings_hist)) < epsilon:
                break

            ## update hist
            np.copyto(p_mappings_hist, p_mappings)

        ## return assignments
        assignments = []
        for row in p_reads.tocsr():
            row = row.toarray().squeeze()
            assignments.append(np.where(row == row.max())[0].tolist())

        ties = defaultdict(set)
        for qseqid, lineage in enumerate(assignments):
            if len(assignment := lineages[lineage]) > 1:
                ties[tuple(assignment)].add(qseqids[qseqid])
            else:
                self.assignments[qseqids[qseqid]] = assignment[0]

        ## resolve ties for equal probability cases using AS, MS and ID
        if ties:
            qset = set.union(*(set(qseqid) for qseqid in ties.values()))
            data = [row for row in data if row[0] in qset]

            for lineages, qseqids in ties.items():
                target = [row for row in data if row[0] in qseqids and row[-1] in lineages]

                scores = defaultdict(lambda: defaultdict(list))
                for row in target:
                    scores[row[-1]]['AS'].append(row[2])
                    scores[row[-1]]['DE'].append(row[3])
                    scores[row[-1]]['ID'].append(row[4])

                ## if all the same, choose the one with known species name
                target = sorted([
                    [
                        np.mean(score['AS']),
                        np.mean(score['DE']),
                        np.mean(score['ID']),
                        not bool(re.search(r' sp\.$| sp\. | sp[0-9]+', lineage.split(';')[-1])),
                        lineage
                    ] for lineage, score in scores.items()
                ], reverse=True)[0][-1]

                for qseqid in qseqids:
                    self.assignments[qseqid] = target

    def run(self, db_kraken=None, skip_profile=False, skip_clean=False,
            max_target_seqs=25, evalue=1e-15, identity=0, subject_cover=75,
            secondary_num=2147483647, secondary_ratio=0.9,
            max_iteration=1000, epsilon=1e-10):
        '''
        Run the pipeline.
        '''
        if db_kraken is not None:
            logger.info('Filtering reads ...')
            self.run_kraken(db_kraken)
            self.parse_kraken()
            logger.info('... removed {} putatively non-prokaryotic reads.'.format(len(self.nset)))

        logger.info('Estimating genome copies ...')
        self.run_diamond(max_target_seqs, evalue, identity, subject_cover)
        self.parse_diamond()
        logger.info('... found {} copies of genomes (bacteria: {}; archaea: {}).'.format(
            sum(self.copies.values()), self.copies['bacteria'], self.copies['archaea']))

        if not skip_profile:
            logger.info('Assigning taxonomy ...')
            self.run_minimap(secondary_num, secondary_ratio)
            self.parse_minimap()
            self.postprocess(max_iteration, epsilon)

            ## fill missing ones according to hits
            replacements = {
                'bacteria': ';'.join(['2|Bacteria'] + ['0|unclassified Bacteria ' + rank for rank in self.ranks[1:]]),
                'archaea': ';'.join(['2157|Archaea'] + ['0|unclassified Archaea ' + rank for rank in self.ranks[1:]])
            }

            ## fit GTDB style
            if self.assignments and '|' not in next(iter(self.assignments.values())).split(';')[0]:
                replacements = {kingdom: ';'.join(
                    rank.split('|')[-1] for rank in replacement.split(';')
                ) for kingdom, replacement in replacements.items()}

            ## count assigned taxonomic labels
            self.hits = [[*hit, self.assignments.get(hit[0], replacements.get(hit[1]))] for hit in self.hits]
            counts, total_counts, lineage2identity = defaultdict(lambda: 0), defaultdict(lambda: 0), defaultdict(list)

            for hit in self.hits:
                total_counts[hit[1]] += 1
                counts[(hit[-1], hit[1])] += 1
                lineage2identity[hit[-1]].append(self.identities.get((hit[0], hit[-1]), (0, 0)))

            copies = defaultdict(lambda: 0)
            for lineage, count in counts.items():
                copies[lineage[0]] += count * self.copies[lineage[1]] / total_counts[lineage[1]]
            total_copies = sum(copies.values())

            # generate a profile output
            self.profile = sorted([
                [*lineage.split(';'), copy, copy / total_copies, *np.mean(lineage2identity.get(lineage), axis=0)] for lineage, copy in copies.items()
            ], key=lambda row: (row[8], row[10], row[7]))

            richness = {'bacteria': 0, 'archaea': 0}
            with open(get_filename(self.file, self.output, '.tsv'), 'w') as w:
                w.write('\t'.join(self.ranks + ['copy', 'abundance', 'identity']) + '\n')

                for row in self.profile:
                    if not re.search('unclassified (Bacteria|Archaea) species', row[6]):
                        richness[row[0].split('|')[-1].lower()] += 1
                    w.write('\t'.join(row[:7] + [f'{row[7]:.3f}', f'{row[8]:e}', f'{row[9]:.4f}/{row[10]:.4f}']) + '\n')

            logger.info('... found {} unique species (bacteria: {}; archaea: {}).'.format(
                sum(richness.values()), richness['bacteria'], richness['archaea']))

        else:
            self.profile = sorted([
                ['Bacteria', self.copies['bacteria'], self.copies['bacteria'] / sum(self.copies.values())],
                ['Archaea', self.copies['archaea'], self.copies['archaea'] / sum(self.copies.values())]
            ], key=lambda row: (row[-2], row[-3]))

            with open(get_filename(self.file, self.output, '.tsv'), 'w') as w:
                w.write('\t'.join(['superkingdom', 'copy', 'abundance']) + '\n')

                for row in self.profile:
                    if row[1] != 0:
                        w.write('\t'.join(row[:1] + [f'{row[1]:.3f}', f'{row[2]:e}']) + '\n')

        ## save reads
        reads = {qseqid: {'remark': 'putatively non-prokaryotic', 'lineage': 'others'} for qseqid in self.nset}
        if not skip_profile:
            reads.update({hit[0]: {
                'remark': 'marker-gene-containing', 'lineage': hit[-1]} for hit in self.hits})
        else:
            reads.update({hit[0]: {
                'remark': 'marker-gene-containing', 'lineage': hit[1].capitalize()} for hit in self.hits})

        with open(get_filename(self.file, self.output, '.json'), 'w') as w:
            json.dump(dict(sorted(reads.items())), w, indent=4)

        ## clean up
        if not skip_clean:
            for file in glob.glob(get_filename(self.file, self.output, '.*.tmp')):
                os.remove(file)
