try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import unittest
from replay_saved_forest import paths_from


class SavedReplayTests(unittest.TestCase):
    def test_paths_and_inverse_transport(self):
        graph={'a':{0:'b'},'b':{3:'a',1:'c'},'c':{4:'b'}}
        paths=paths_from(graph,'a')
        self.assertEqual(paths,{'a':[],'b':[0],'c':[0,1]})
        transport=[(m+3)%6 for m in reversed(paths['c'])]
        x='c'
        for m in transport:x=graph[x][m]
        self.assertEqual(x,'a')

    def test_disjoint_partial_graph_is_not_claimed_same_component(self):
        graph={'a':{0:'b'},'b':{3:'a'},'c':{1:'d'},'d':{4:'c'}}
        self.assertNotIn('c',paths_from(graph,'a'))


if __name__=='__main__':unittest.main()
