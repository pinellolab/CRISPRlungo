import os

os.system('rm testAgreement.settings.CRISPRlungo.*')
stream = os.popen('python ../../CRISPRlungo.py testAgreement.settings')
output = stream.read()
print(output)
test_count = 0

test_col = 14

with open("testAgreement.settings.CRISPRlungo.final.final_assignments",'r') as fin:
    head = fin.readline()
    head_els = head.split("\t")
    assert head_els[test_col] == 'r1_r2_support_status', "Expecting 'r1_r2_support_status' in column " + str(test_col) + "(got " + head_els[test_col] + ")"
    for line in fin:
        line_els = line.split("\t")
        true_status_els = line_els[0].split("#")
        true_status = true_status_els[1]
        computed_status = line_els[test_col].replace(' ','_')
        line_id = true_status_els[0]
        print("line id: " + line_id)
        print('true status: ' + true_status)
        print('computed status: ' + computed_status)

        assert true_status == computed_status, "For " + line_id + " got status " + computed_status + " (Expecting " + true_status + ")"
        test_count += 1

print(str(test_count) + " tests passed")