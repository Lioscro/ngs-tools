import os
from unittest import TestCase

from ngs_tools import gtf

from . import mixins


class TestGtf(mixins.TestMixin, TestCase):

    def test_read(self):
        g = gtf.Gtf(self.gtf_path, 'r')
        entries = list(g)
        self.assertEqual(9, len(entries))
        self.assertEqual('2', entries[0].chromosome)
        self.assertEqual('exon', entries[0].feature)
        self.assertEqual(2, entries[0].start)
        self.assertEqual(3, entries[0].end)
        self.assertEqual('-', entries[0].strand)
        self.assertEqual({
            'gene_id': 'GENE_C',
            'gene_name': 'GENE_C_NAME',
            'transcript_id': 'TRANSCRIPT_C',
            'gene_source': 'havana',
            'gene_biotype': 'TEC',
        }, entries[0].attributes)

    def test_write(self):
        path = os.path.join(self.temp_dir, 'test.gtf')
        g = gtf.Gtf(path, 'w')
        line = (
            '2\thavana\texon\t2\t3\t.\t-\t.\t'
            'gene_id "GENE_C"; gene_name "GENE_C_NAME"; '
            'transcript_id "TRANSCRIPT_C"; gene_source "havana"; gene_biotype "TEC";'
        )
        g.write(gtf.GtfEntry(line))
        g.close()
        with open(path, 'r') as f:
            self.assertEqual(f'{line}\n', f.read())

    def test_parse_gtf(self):
        entries = list(gtf.parse_gtf(self.gtf_path))
        self.assertEqual(9, len(entries))
        self.assertEqual('2', entries[1].chromosome)
        self.assertEqual('exon', entries[1].feature)
        self.assertEqual(10, entries[1].start)
        self.assertEqual(14, entries[1].end)
        self.assertEqual('-', entries[1].strand)
        self.assertEqual({
            'gene_id': 'GENE_C',
            'gene_name': 'GENE_C_NAME',
            'transcript_id': 'TRANSCRIPT_C',
            'gene_source': 'havana',
            'gene_biotype': 'TEC',
        }, entries[1].attributes)

    def test_parse_gtf_with_filter_func(self):
        entries = list(
            gtf.parse_gtf(
                self.gtf_path,
                filter_func=lambda entry: entry.feature == 'exon'
            )
        )
        self.assertEqual(5, len(entries))

    def test_genes_and_transcripts_from_gtf(self):
        gene_infos, transcript_infos = gtf.genes_and_transcripts_from_gtf(
            self.gtf_path
        )
        self.assertEqual({
            'GENE_D': {
                'transcripts': ['GENE_D'],
                'segment': gtf.Segment(0, 2),
                'chromosome': '1',
                'strand': '+',
                'gene_name': 'GENE_D_NAME'
            },
            'GENE_C': {
                'transcripts': ['TRANSCRIPT_C'],
                'segment': gtf.Segment(1, 14),
                'chromosome': '2',
                'strand': '-',
                'gene_name': 'GENE_C_NAME'
            },
            'GENE_B': {
                'transcripts': ['TRANSCRIPT_B'],
                'segment': gtf.Segment(0, 10),
                'chromosome': '2',
                'strand': '+',
                'gene_name': 'GENE_B_NAME'
            },
            'GENE_A': {
                'transcripts': ['TRANSCRIPT_A'],
                'segment': gtf.Segment(0, 2),
                'chromosome': '1',
                'strand': '+',
                'gene_name': 'GENE_A_NAME'
            }
        }, gene_infos)
        self.assertEqual({
            'GENE_D': {
                'gene_id': 'GENE_D',
                'transcript_name': '',
                'segment': gtf.Segment(0, 2),
                'exons': gtf.SegmentCollection([gtf.Segment(0, 2)]),
                'introns': gtf.SegmentCollection()
            },
            'TRANSCRIPT_C': {
                'gene_id':
                    'GENE_C',
                'transcript_name':
                    '4933401J01Rik-201',
                'segment':
                    gtf.Segment(1, 14),
                'exons':
                    gtf.SegmentCollection([
                        gtf.Segment(1, 3),
                        gtf.Segment(9, 14)
                    ]),
                'introns':
                    gtf.SegmentCollection([gtf.Segment(3, 9)])
            },
            'TRANSCRIPT_B': {
                'gene_id':
                    'GENE_B',
                'transcript_name':
                    '4933401J01Rik-201',
                'segment':
                    gtf.Segment(0, 10),
                'exons':
                    gtf.SegmentCollection([
                        gtf.Segment(1, 2), gtf.Segment(4, 8)
                    ]),
                'introns':
                    gtf.SegmentCollection([
                        gtf.Segment(0, 1),
                        gtf.Segment(2, 4),
                        gtf.Segment(8, 10)
                    ])
            },
            'TRANSCRIPT_A': {
                'gene_id': 'GENE_A',
                'transcript_name': '4933401J01Rik-201',
                'segment': gtf.Segment(0, 2),
                'exons': gtf.SegmentCollection([gtf.Segment(0, 2)]),
                'introns': gtf.SegmentCollection(),
            }
        }, transcript_infos)

    def test_genes_and_transcripts_from_gtf_zero_length(self):
        gene_infos, transcript_infos = gtf.genes_and_transcripts_from_gtf(
            self.zero_length_gtf_path
        )
        self.assertEqual({
            'Arnt': {
                'segment': gtf.Segment(183054636, 183056584),
                'chromosome': 'NC_051337.1',
                'strand': '+',
                'gene_name': '',
                'transcripts': ['NM_012780.2']
            }
        }, gene_infos)
        self.assertEqual({
            'NM_012780.2': {
                'gene_id':
                    'Arnt',
                'segment':
                    gtf.Segment(183054636, 183056584),
                'transcript_name':
                    '',
                'exons':
                    gtf.SegmentCollection([
                        gtf.Segment(183054636, 183055283),
                        gtf.Segment(183055283, 183056584)
                    ]),
                'introns':
                    gtf.SegmentCollection()
            }
        }, transcript_infos)


class TestSegment(mixins.TestMixin, TestCase):

    def test_init(self):
        segment = gtf.Segment(3, 4)
        self.assertEqual(3, segment.start)
        self.assertEqual(4, segment.end)
        self.assertEqual(1, segment.width)
        with self.assertRaises(gtf.SegmentError):
            gtf.Segment(1, 0)

    def test_is_in(self):
        segment = gtf.Segment(0, 10)
        self.assertTrue(segment.is_in(0))
        self.assertFalse(segment.is_in(10))

    def test_is_exclusive(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 20)
        self.assertTrue(segment1.is_exclusive(segment3))
        self.assertFalse(segment1.is_exclusive(segment2))

    def test_is_overlapping(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 20)
        self.assertFalse(segment1.is_overlapping(segment3))
        self.assertTrue(segment1.is_overlapping(segment2))

    def test_is_subset(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(0, 5)
        segment4 = gtf.Segment(5, 15)
        self.assertTrue(segment2.is_subset(segment1))
        self.assertTrue(segment3.is_subset(segment1))
        self.assertFalse(segment4.is_subset(segment1))

    def test_is_superset(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(0, 5)
        segment4 = gtf.Segment(5, 15)
        self.assertTrue(segment1.is_superset(segment2))
        self.assertTrue(segment1.is_superset(segment3))
        self.assertFalse(segment1.is_superset(segment4))

    def test_comparison(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(0, 10)
        segment3 = gtf.Segment(5, 10)
        segment4 = gtf.Segment(0, 5)
        self.assertTrue(segment1 == segment2)
        self.assertTrue(segment1 < segment3)
        self.assertTrue(segment1 > segment4)

    def test_flank(self):
        segment = gtf.Segment(0, 10)
        flanked = segment.flank(2)
        self.assertEqual(0, flanked.start)
        self.assertEqual(12, flanked.end)

        flanked = segment.flank(2, left=1, right=10)
        self.assertEqual(1, flanked.start)
        self.assertEqual(10, flanked.end)


class TestSegmentCollection(mixins.TestMixin, TestCase):

    def test_init(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(5, 15)
        collection = gtf.SegmentCollection(segments=[segment1, segment2])
        self.assertEqual(0, collection.start)
        self.assertEqual(15, collection.end)
        self.assertEqual(1, len(collection))
        self.assertEqual(gtf.Segment(0, 15), collection.segments[0])

    def test_bool(self):
        collection1 = gtf.SegmentCollection()
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(0, 10)])
        self.assertFalse(bool(collection1))
        self.assertTrue(bool(collection2))

    def test_add_segment(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 15)
        collection = gtf.SegmentCollection(segments=[segment1, segment3])
        collection.add_segment(segment2)
        self.assertEqual(3, len(collection))
        self.assertEqual([segment1, segment2, segment3], collection.segments)

    def test_add_collection(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(10, 15)
        collection1 = gtf.SegmentCollection(segments=[segment1, segment3])
        collection2 = gtf.SegmentCollection(segments=[segment2])
        collection1.add_collection(collection2)
        self.assertEqual(3, len(collection1))
        self.assertEqual([segment1, segment2, segment3], collection1.segments)

    def test_invert(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(10, 15)
        collection = gtf.SegmentCollection(segments=[segment1, segment2])
        with self.assertRaises(gtf.SegmentCollectionError):
            collection.invert(gtf.Segment(10, 30))

        self.assertEqual(
            gtf.SegmentCollection([gtf.Segment(5, 10),
                                   gtf.Segment(15, 20)]),
            collection.invert(gtf.Segment(0, 20))
        )

    def test_collapse(self):
        collection = gtf.SegmentCollection()
        collection._segments = [gtf.Segment(0, 10), gtf.Segment(5, 15)]
        collection.collapse()
        self.assertEqual(1, len(collection))
        self.assertEqual(gtf.Segment(0, 15), collection.segments[0])

    def test_span_is_exclusive(self):
        collection1 = gtf.SegmentCollection(
            segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)]
        )
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(5, 10)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(15, 20)])
        self.assertFalse(collection1.span_is_exclusive(collection2))
        self.assertTrue(collection1.span_is_exclusive(collection3))

    def test_is_overlapping(self):
        collection1 = gtf.SegmentCollection(
            segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)]
        )
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(5, 10)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(0, 10)])
        self.assertFalse(collection1.is_overlapping(collection2))
        self.assertTrue(collection1.is_overlapping(collection3))

    def test_is_subset(self):
        collection1 = gtf.SegmentCollection(
            segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)]
        )
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(3, 13)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(1, 2)])
        self.assertFalse(collection2.is_subset(collection1))
        self.assertTrue(collection3.is_subset(collection1))

    def test_is_superset(self):
        collection1 = gtf.SegmentCollection(
            segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)]
        )
        collection2 = gtf.SegmentCollection(segments=[gtf.Segment(3, 13)])
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(1, 2)])
        self.assertFalse(collection1.is_superset(collection2))
        self.assertTrue(collection1.is_superset(collection3))

    def test_from_positions(self):
        positions = [0, 1, 2, 3, 3, 5, 6, 7, 8]
        collection = gtf.SegmentCollection.from_positions(positions)
        self.assertEqual([gtf.Segment(0, 4),
                          gtf.Segment(5, 9)], collection.segments)

    def test_from_collections(self):
        segment1 = gtf.Segment(0, 5)
        segment2 = gtf.Segment(5, 10)
        segment3 = gtf.Segment(9, 15)
        collection1 = gtf.SegmentCollection(segments=[segment1, segment3])
        collection2 = gtf.SegmentCollection(segments=[segment2])
        collection = gtf.SegmentCollection.from_collections(
            collection1, collection2
        )
        self.assertEqual([gtf.Segment(0, 5),
                          gtf.Segment(5, 15)], collection.segments)

    def test_comparison(self):
        collection1 = gtf.SegmentCollection(
            segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)]
        )
        collection2 = gtf.SegmentCollection(
            segments=[gtf.Segment(0, 5), gtf.Segment(10, 15)]
        )
        collection3 = gtf.SegmentCollection(segments=[gtf.Segment(1, 2)])
        self.assertTrue(collection1 == collection2)
        self.assertFalse(collection1 == collection3)

    def test_flank(self):
        segment1 = gtf.Segment(0, 10)
        segment2 = gtf.Segment(12, 15)
        segment3 = gtf.Segment(20, 25)
        collection = gtf.SegmentCollection(
            segments=[segment1, segment2, segment3]
        )
        flanked = collection.flank(2)
        self.assertEqual(2, len(flanked))
        self.assertEqual(0, flanked[0].start)
        self.assertEqual(17, flanked[0].end)
        self.assertEqual(18, flanked[1].start)
        self.assertEqual(27, flanked[1].end)
