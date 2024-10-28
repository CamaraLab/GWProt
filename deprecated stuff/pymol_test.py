import pymol

def test_run():
    pymol.finish_launching(["pymol", "-pc"])
    print('launched')
    pymol.cmd.fetch('1hq0')
    print('fetched')
    pymol.cmd.transform_selection('1hq0', [1,0,0,5,0,1,0,5,0,0,1,5,0,0,0,0])
    print('transformed, done!')

