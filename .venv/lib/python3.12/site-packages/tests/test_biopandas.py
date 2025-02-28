def test_featurizer():
    from biopandas.pdb import PandasPdb

    file_unprotonated = "data/pdb/rcsb/1ahw.pdb"
    ppdb = PandasPdb().read_pdb(file_unprotonated)

    print(ppdb.df['ATOM'].head().to_string())

test_featurizer()