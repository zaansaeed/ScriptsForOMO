def test_occupancy_claudio():
    from protkit.file_io import PDBIO
    prot = PDBIO.load("data/claudio/3s9d.pdb")[0]
    prot.fix_disordered_atoms()

    print(prot)

test_occupancy_claudio()