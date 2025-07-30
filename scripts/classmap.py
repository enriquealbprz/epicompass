# Import libraries
import csv
import sys
import os
import argparse

class Sample: # Object 'sample' and its information
    def __init__(self,group,bssid,info,name,age,sex,lifestage):
        self.bssid = bssid
        self.name = name
        self.group = group
        self.info = info
        self.age = age
        self.sex = sex
        self.lifestage = lifestage

    def __str__(self):
        return f"[{self.bssid}] {self.name} ({self.group}, {self.info})"

# Read metadata
SD = open("../data/main_metadata_table.tsv")
tsv_reader = csv.reader(SD, delimiter="\t")
next(tsv_reader)  # Skip header
sample_dict = {}
for line in tsv_reader:
    (id,g,ct,secondary,i,ori,pert,life,ag,ageunits,sx,type,col,cat,proj,don,n) = line
    sample = Sample(f"{g}", f"{id}", f"{i}", f"{n}", f"{ag}", f"{sx}", f"{life}")
    sample_dict.update({f"{id}": sample})

# Parse arguments
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("criterion", help="sample classification criterion ('tissue', 'age', 'sex' or 'lifestage')", type=str)
    parser.add_argument("--output", help="Output file for classmap (.tsv)", type=str, default='../data/classmap.tsv')    

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Generate output file according to the chosen criterion
    if not args.output.endswith('.tsv'):
        raise ValueError("Output file extension must be .tsv")

    if args.criterion == "tissue":
        with open(args.output, "w") as f:
            for s in sample_dict:
                f.write(f"{sample_dict[s].bssid}\t{sample_dict[s].group}\n")
    elif args.criterion == "age":
        with open(args.output, "w") as f:
            for s in sample_dict:
                f.write(f"{sample_dict[s].bssid}\t{sample_dict[s].age}\n")
    elif args.criterion == "sex":
        with open(args.output, "w") as f:
            for s in sample_dict:
                f.write(f"{sample_dict[s].bssid}\t{sample_dict[s].sex}\n")
    elif args.criterion == "lifestage":
        with open(args.output, "w") as f:
            for s in sample_dict:
                f.write(f"{sample_dict[s].bssid}\t{sample_dict[s].lifestage}\n")
    elif args.criterion != "tissue" and args.criterion != "age" and args.criterion != "sex" and args.criterion != "lifestage":
        raise ValueError("Unexpected criterion; must be 'tissue', 'age', 'sex' or 'lifestage'.")

if __name__ == "__main__":
    main()