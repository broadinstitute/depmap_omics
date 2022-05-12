import re
import sys
import bgzip
import gzip
import collections

print(sys.argv)

line_buff = collections.deque([])
line_pos = collections.deque([])
check_index = 0
size_tocheck = int(sys.argv[3])

with gzip.open(sys.argv[1], "r+") as f:
    with open(sys.argv[2], "wb") as raw_out:
        with bgzip.BGZipWriter(raw_out) as f_out:
            for i, line in enumerate(f):
                # print progress
                if i % 100_000 == 0:
                    print("Processed " + str(i) + " variants")
                original_string = line.decode("utf-8")
                # process header
                if original_string[0] == "#":
                    f_out.write(original_string.encode())
                else:
                    new_string = re.sub(
                        r"AS_FilterStatus=(.*?);",
                        lambda x: "AS_FilterStatus="
                        + x.group(1)
                        .replace("|", "~")
                        .replace(",", "|")
                        .replace("~", ",")
                        + ";",
                        original_string,
                    )
                    to_print = new_string.split("\t")
                    to_print[7] = to_print[7].replace(" ", "")

                    # Add to line buffer
                    line_buff.append(to_print)
                    current_pos = int(to_print[1])
                    line_pos.append(current_pos)

                    # Check lines > size_tocheck from current position
                    while abs(line_pos[check_index] - current_pos) > size_tocheck:
                        if line_buff[check_index][6] in [
                            "clustered_events",
                            "clustered_events;haplotype",
                        ]:
                            n_somatic = 0
                            for neighbor in line_buff:
                                if "germline" not in neighbor[6] and abs(
                                    int(neighbor[1]) - line_pos[check_index]
                                ) <= int(size_tocheck / 2):
                                    n_somatic += 1
                            if n_somatic <= 2:
                                line_buff[check_index][6] = "PASS"
                        check_index += 1

                    # Write lines >2*size_tocheck from current position
                    while abs(line_pos[0] - current_pos) > 2 * size_tocheck:
                        line_to_write = line_buff.popleft()
                        to_print = "\t".join(line_to_write)
                        f_out.write(to_print.encode())

                        # roll back the queue
                        check_index -= 1
                        line_pos.popleft()
            # Write the final lines
            for line_to_write in line_buff:
                to_print = "\t".join(line_to_write)
                f_out.write(to_print.encode())
