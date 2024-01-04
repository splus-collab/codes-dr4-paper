#!/bin/usr/python3

import json
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create authors.tex file from list_authors.json')
    parser.add_argument('--workdir', type=str, default=os.getcwd(),
                        help='Working directory. Default: current directory')
    parser.add_argument('--authors', type=str, default='config/list_authors.json',
                        help='Authors file. Default: config/list_authors.json')

    args = parser.parse_args()
    return args


def get_authors_list(args):
    """Get authors list from list_authors.json file"""
    authors_json = json.load(
        open(os.path.join(args.workdir, args.authors), 'r'))

    authors = authors_json['authors']

    # authors = sorted(authors.items(), key=lambda x: (
    #     x[1].get('position', 99), x[1]['fullname']))
    authors = sorted(authors.items(), key=lambda x: (
        x[1]['position'], x[1]['fullname']))

    return authors_json, authors


def prepare_authors_files(args, authors_json, authors):
    """Prepare authors.txt and authors.tex files"""
    institutions = []
    for author in authors:
        for institution in author[1]['institutions']:
            if institution not in institutions:
                institutions.append(institution)

    a = open(os.path.join(args.workdir,
             f'extras/{args.authors.split("/")[-1].strip(".json")}.txt'), 'w')
    with open(os.path.join(args.workdir, f'extras/{args.authors.split("/")[-1].strip(".json")}.tex'), 'w') as f:
        f.write('\\author{')
        i = 1
        listemails = []
        for author in authors:
            a.write('%i - %s\n' % (i, author[1]['fullname']))
            listemails.append(f'<{author[1]["email"]}>')
            i += 1
            if author is authors[-1]:
                endofline = '\n'
            else:
                endofline = ',\n'
            insts = [institutions.index(inst) + 1
                     for inst in author[1]['institutions']]
            if author[1]['orcid'] != "":
                f.write('%s\\inst{%s}\\orcid{%s}%s' % (
                    author[0], ','.join(str(x) for x in sorted(insts)), author[1]['orcid'], endofline))
            else:
                f.write('%s\\inst{%s}%s' %
                        (author[0], ','.join(str(x) for x in sorted(insts)), endofline))
        f.write('}\n')
        f.write('\\institute{')
        for institution in institutions:
            if institution is institutions[0]:
                endofline = '\\\\\\email{herpich@ast.cam.ac.uk}\\and\n'
            elif institution is institutions[-1]:
                endofline = '\n}'
            else:
                endofline = '\\and\n'

            f.write('%s, %s, %s, %s, %s%s' %
                    (authors_json['institutions'][institution]['name'],
                     authors_json['institutions'][institution]['addressline'],
                     authors_json['institutions'][institution]['city'],
                     authors_json['institutions'][institution]['postcode'],
                     authors_json['institutions'][institution]['country'],
                     endofline))

        f.write('\n')
        f.write('\n')

        f.write('\\section*{aknowledgements}\n')
        for author in authors:
            if author[1]['acknowledgements'] != "":
                f.write('%s\n' % author[1]['acknowledgements'])

        f.write('\n')
        f.write('\n')

        f.write('\\section{Authors Affiliations}\\label{ap:affiliations}\n')
        for institution in institutions:
            f.write('\n')
            f.write('\\noindent\n')
            f.write('$^{%i}$ %s, %s, %s, %s, %s\n' %
                    (institutions.index(institution) + 1,
                     authors_json['institutions'][institution]['name'],
                     authors_json['institutions'][institution]['addressline'],
                     authors_json['institutions'][institution]['city'],
                     authors_json['institutions'][institution]['postcode'],
                     authors_json['institutions'][institution]['country']))

    a.write('\n\n')
    for author_email in listemails:
        a.write(f'{author_email}, ')
    a.close()


def main():
    args = parse_args()

    authors_json, authors = get_authors_list(args)

    prepare_authors_files(args, authors_json, authors)

    return


if __name__ == '__main__':
    main()
