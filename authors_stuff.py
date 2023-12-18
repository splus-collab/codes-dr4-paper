#!/bin/usr/python3

import json

authors_json = json.load(open('config/list_authors.json', 'r'))

authors = authors_json['authors']
# order authors by the key position knowing that authors is a dict
authors = sorted(authors.items(), key=lambda x: x[1]['position'])
# each author has a key 'institution' that is a list of institutions. Need to order them by apperance and enumerate them starting in 1
institutions = []
for author in authors:
    for institution in author[1]['institutions']:
        if institution not in institutions:
            institutions.append(institution)

a = open('authors_names.txt', 'w')
with open('authors.tex', 'w') as f:
    f.write('\\author{')
    for author in authors:
        a.write('%s\n' % author[1]['fullname'])
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


a.close()
