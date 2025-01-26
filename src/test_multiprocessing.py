from multiprocessing import Pool

def read_pairs():
    with open("./data/in.txt", "r") as f:
        while True:
            l1 = f.readline()
            l2 = f.readline()
            if not l1 or not l2:
                break
            yield l1, l2


def process_pair(pair):
    l1, l2 = pair
    with open("./data/out.txt", "a") as out:
        out.write(l1.strip() + " " + l2.strip() + "\n")

if __name__ == "__main__":
    with open("./data/out.txt", "w") as f:
        f.write("Output\n")

    with Pool(8) as pool:
        pool.map(process_pair, read_pairs())
