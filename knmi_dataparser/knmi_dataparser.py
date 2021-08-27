import csv

# date of simulation in YYYYMMDD
simulation_date = "20190810"


def simple_forcing(file_name):

    with open(file_name, 'r') as fh:
        line = " "
        while line[0] != '#':
            line = fh.readline()

        date_idx = 1
        hh_idx = 2
        t_idx = 7
        u_idx = 17

        header = line.split(',')
        header = [c.strip() for c in header]

        assert(header[date_idx] == "YYYYMMDD")
        assert(header[hh_idx] == "HH")
        assert(header[t_idx] == "T")
        assert(header[u_idx] == "U")

        fh.readline()

        data = csv.reader(fh, delimiter=',')

        print("Hourly temperature and humidity for {}".format(simulation_date))
        print("Var/Ti\tT\t\tq")

        for row in data:
            if row[date_idx] == simulation_date:
                print("{}:00\t{}\t{}".format(int(row[hh_idx])-1, int(row[t_idx])/10.0, row[u_idx]))


def main():
    simple_forcing("input/uurgeg_269_2011-2020.txt")


if __name__ == '__main__':
    main()
