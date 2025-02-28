def test_reduce():
    from protkit.file_io.pdb_io import PDBIO
    from protkit.tools.reduce_adaptor import ReduceAdaptor

    file_unprotonated = "data/pdb/rcsb/1ahw.pdb"
    file_temp_protonated = "data/pdb/rcsb/1ahw_reduce_h.pdb"
    file_protonated = "data/pdb/rcsb/1ahw_h.pdb"

    protein = PDBIO.load(file_unprotonated)[0]
    reduce = ReduceAdaptor("/usr/local/bin/reduce", quiet=True)

    protein2 = reduce.protonate(protein, reduce_output_pdb_file_path=file_temp_protonated)
    print(protein2)
    PDBIO.save(protein2, file_protonated)


test_reduce()
