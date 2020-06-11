import argparse


def check(val, expected):
    rel_tolerance = 10e-5
    abs_tolerance = 10e-5
    return abs(val - expected) <= max(rel_tolerance * max(abs(val), abs(expected)), abs_tolerance)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile', help='Output file produced by your program')
    parser.add_argument("correctfile", help="File containing correct answers")
    args = parser.parse_args()

    with open(args.outfile, 'r') as outfile:
        with open(args.correctfile, 'r') as correctfile:
            out_lines = list(filter(lambda l : l != "\n", outfile))
            correct_lines = list(filter(lambda l : l != "\n", correctfile))

            if len(out_lines) != len(correct_lines):
                print("[%s] different number of lines" % (args.correctfile))
                exit(1)
            for i in range(len(out_lines)):
                out_line = out_lines[i].split(' ')
                correct_line = correct_lines[i].split(' ')
                for j in range(len(out_line)):
                    if not check(float(out_line[j]), float(correct_line[j])):
                        print("[%s] line %d column %d expected %s, actual %s" %
                              (args.correctfile, i, j, correct_line[j], out_line[j]))
                        exit(1)
