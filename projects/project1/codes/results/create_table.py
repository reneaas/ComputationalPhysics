from os import system

def main():
    tableCells = int(input("Enter the number of tablecells\n\t"))
    dataFile = input("Enter the name of the datafile\n\t")
    data = read_data(dataFile)
    print(data)
    dataLength = len(data)
    with open("table.tex", "w") as f:
        f.write("\\begin{tabular}{")
        for i in range(0, tableCells):
            if i + 1 == tableCells:
                f.write("|c|")
            else:
                f.write("|c")
        f.write("}\n")
        for x in range(0, dataLength):
            f.write("\t\hline\n")
            data[x] = data[x].replace(" ", "&")
            f.write("\t" + data[x] + "\\\\\n")
        f.write("\hline\n")
        f.write("\\end{tabular}")


def read_data(dataFile):
    with open(dataFile, "r") as df:
            data = df.read()

    data = data.split("\n")
    return data


main()
