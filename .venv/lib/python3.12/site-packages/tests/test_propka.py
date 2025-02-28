def test_propka():
    from protkit.file_io.pdb_io import PDBIO
    from protkit.file_io.pqr_io import PQRIO
    from protkit.tools.propka_adaptor import PropkaAdaptor

    # file_unprotonated = "data/pdb/rcsb/1ahw.pdb"
    # file_unprotonated = "data/pdb/rcsb/1a4y.pdb"
    # file_unprotonated = "data/pdb/rcsb/3i40.pdb"
    # file_unprotonated = "data/pdb/rcsb/4nkq.pdb"
    file_unprotonated = "data/pdb/rcsb/6bom.pdb"

    protein = PDBIO.load(file_unprotonated)[0]
    propka = PropkaAdaptor(ph=7.2)

    protein2 = propka.calculate_pka(protein)
    print(protein2)


test_propka()