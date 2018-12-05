from prakt.fasta_parser.fasta_parser import parse_fasta

def test_parse_fast():
    records = parse_fasta("data/test.fasta")

    assert len(records) == 2
    
    record1 = records[0]
    record2 = records[1]

    assert str(record1.seq) == "FancySequenceA"
    assert str(record2.seq) == "FancysequenceB"

    assert record1.id == "idA"
    assert record2.id == "idB"