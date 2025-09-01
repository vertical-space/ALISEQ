import simulation

def test_fasta():
    obs = simulation.fasta({'id1':'atcg'})
    exp = """>id1
atcg"""
    assert obs == exp


def newFunction():
    '''I created this new fuction to test some changes...
    '''
    print("testing my changes")
