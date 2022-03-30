import os

os.system('rm testCustomTargetPos.settings.CRISPRlungo.*')
stream = os.popen('python ../../CRISPRlungo.py testCustomTargetPos.settings')
output = stream.read()
print(output)
test_count = 0

with open("testCustomTargetPos.settings.CRISPRlungo.final.final_assignments",'r') as fin:
    head = fin.readline()
    for line in fin:
        line_els = line.split("\t")
        true_loc_els = line_els[0].split("_")
        aligned_locs = line_els[4].split("~")
        line_id = true_loc_els.pop(0)
        anno = true_loc_els.pop(0)
        print("line id: " + line_id)
        print('true locs: ' + str(true_loc_els))
        print('aligned losc: ' + str(aligned_locs))
        if len(true_loc_els) < 4:
            true_chr,true_start,true_end = true_loc_els
            assert len(aligned_locs) == 1, "For " + line_id + " aligned to more than one loc (true " + str(true_loc_els) + " vs " + str(aligned_locs) + ")"
            aligned_chr,aligned_pos = aligned_locs[0].split(":")
            aligned_start,aligned_end = aligned_pos.split("-")

            assert true_chr == aligned_chr, "For " + line_id + " got mismatched chr ( true " + true_chr + " vs " + aligned_chr + ")"
            assert true_start == aligned_start, "For " + line_id + " got mismatched start ( true " + true_start + " vs " + aligned_start + ")"
            assert true_end == aligned_end, "For " + line_id + " got mismatched end ( true " + true_end + " vs " + aligned_end + ")"

        else:
            assert len(aligned_locs) == 2, "For " + line_id + " aligned to more than one loc (true " + str(true_loc_els) + " vs " + str(aligned_locs) + ")"

            true_chr1,true_start1,true_end1, true_chr2,true_start2,true_end2 = true_loc_els

            aligned_chr1,aligned_pos1 = aligned_locs[0].split(":")
            aligned_start1,aligned_end1 = aligned_pos1.split("-")
            aligned_chr2,aligned_pos2 = aligned_locs[1].split(":")
            aligned_start2,aligned_end2 = aligned_pos2.split("-")

            assert true_chr1 == aligned_chr1, "For " + line_id + " got mismatched chr1 ( true " + true_chr1 + " vs " + aligned_chr1 + ")"
            assert true_start1 == aligned_start1, "For " + line_id + " got mismatched start1 ( true " + true_start1 + " vs " + aligned_start1 + ")"
            assert true_end1 == aligned_end1, "For " + line_id + " got mismatched end1 ( true " + true_end1 + " vs " + aligned_end1 + ")"
            assert true_chr2 == aligned_chr2, "For " + line_id + " got mismatched chr2 ( true " + true_chr2 + " vs " + aligned_chr2 + ")"
            assert true_start2 == aligned_start2, "For " + line_id + " got mismatched start2 ( true " + true_start2 + " vs " + aligned_start2 + ")"
            assert true_end2 == aligned_end2, "For " + line_id + " got mismatched end2 ( true " + true_end2 + " vs " + aligned_end2 + ")"
        test_count += 1

print(str(test_count) + " tests passed")
