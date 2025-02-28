def test_pdb2pqr():
    from protkit.file_io.pdb_io import PDBIO
    from protkit.file_io.pqr_io import PQRIO
    from protkit.tools.pdb2pqr_adaptor import PDB2PQRAdaptor

    file_unprotonated = "data/pdb/rcsb/1ahw.pdb"
    file_pdb2pqr_pdb = "data/pdb/rcsb/1ahw_pdb2pqr.pdb"
    file_pdb2pqr_pqr = "data/pqr/rcsb/1ahw_pdb2pqr.pqr"
    file_pqr = "data/pqr/rcsb/1ahw.pqr"

    protein = PDBIO.load(file_unprotonated)[0]
    pdb2pqr = PDB2PQRAdaptor()

    protein2 = pdb2pqr.run(protein, output_pqr_file_path=file_pdb2pqr_pqr, output_pdb_file_path=file_pdb2pqr_pdb)
    print(protein2)
    PQRIO.save(protein2, file_pqr)


test_pdb2pqr()