import argparse
import trios

argparser = argparse.ArgumentParser(description = 'Pulls trios out of the PC-relate/KING output.')
argparser.add_argument('--kinship', metavar = 'file', dest = 'in_kinship_file', required = True, help = 'File with estimated kinship coefficients (KING or PC-relate format).')
argparser.add_argument('--format', metavar = 'name', dest = 'in_format', required = True, choices = ['pc-relate', 'king'], help = 'File format: pc-relate, king')
argparser.add_argument('--sex', metavar = 'file', dest = 'in_sex_file', required = True, help = 'File with two tab-separated columns: id, sex (M/F).')
argparser.add_argument('--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')

if __name__ == '__main__':
   args = argparser.parse_args()

   if args.in_format == 'pc-relate':
      kinship_entries = trios.pcrelate_reader(args.in_kinship_file)
   elif args.in_format == 'king':
      kinship_entries = trios.king_reader(args.in_kinship_file)
   else:
       raise Exception('Specified format %s is not supported.' % args.in_format)

   trios.write_trios(trios.get_trios(kinship_entries, trios.sex_reader(args.in_sex_file), 3), args.out_file)


