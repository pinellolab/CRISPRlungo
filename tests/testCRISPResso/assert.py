import os

os.system('rm -rf testCRISPResso.settings.CRISPRlungo.*')
stream = os.popen('python ../../CRISPRlungo.py testCRISPResso.settings')
output = stream.read()
print(output)
test_count = 0

test_classification_col = 2
test_crispresso_col = 10

with open("testCRISPResso.settings.CRISPRlungo.CRISPResso.annotated_final_assignments",'r') as fin:
    head = fin.readline()
    head_els = head.split("\t")
    assert head_els[test_classification_col] == 'classification', "Expecting 'classification' in column " + str(test_classification_col) + "(got " + head_els[test_classification_col] + ")"
    for line in fin:
        line_els = line.strip().split("\t")
        expected_status_els = line_els[0].split("_")
        expected_classification = expected_status_els[1]
        expected_crispresso = " ".join(expected_status_els[2:4])
        computed_classification = line_els[test_classification_col]
        computed_crispresso = line_els[test_crispresso_col]

        line_id = expected_status_els[0]
        print("line id: " + line_id)
        print('expected classification status: ' + expected_classification)
        print('computed classification status: ' + computed_classification)


        assert expected_classification == computed_classification, "For " + line_id + " got classification " + computed_classification + " (Expecting " + expected_classification + ")"

        print('expected crispresso status: ' + expected_crispresso)
        print('computed crispresso status: ' + computed_crispresso)
        assert expected_crispresso == computed_crispresso, "For " + line_id + " got crispresso " + computed_crispresso + " (Expecting " + expected_crispresso + ")"
        test_count += 1

print(str(test_count) + " tests passed")
