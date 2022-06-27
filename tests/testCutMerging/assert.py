import os

os.system('rm -rf testCutMerging.settings.CRISPRlungo.*')
stream = os.popen('python ../../CRISPRlungo.py testCutMerging.settings')
output = stream.read()
print(output)
test_count = 0

with open("testCutMerging.settings.CRISPRlungo.final.final_cut_points.txt",'r') as fin:
    head = fin.readline()
    head_els = head.split("\t")

    line1 = fin.readline()
    assert line1.strip() == 'chr1\t900049\t2\tNovel'
    line2 = fin.readline()
    assert line2.strip() == 'chr1\t900183\t2\tNovel'
    line3 = fin.readline()
    assert line3.strip() == 'chr5\t500050\t1\tOn-target'


print(str(test_count) + " tests passed")
