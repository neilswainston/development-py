'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.utils import plate_utils


def _write_overview(designs, filename):
    '''Writes an overview file.'''
    overview_file = open(filename + '_overview.xls', 'w+')

    domino_pairs = {}

    line = '\t'.join(['Plasmid id',
                      'Part 1 id',
                      'Part 2 id',
                      'Domino id',
                      'Domino seq',
                      'Domino 1 seq',
                      'Domino 1 Tm',
                      'Domino 2 seq',
                      'Domino 2 Tm'])

    overview_file.write(line + '\n')

    for design_id, design in designs.iteritems():
        for domino in design['dominoes']:
            pair = list(domino[0])
            seq = domino[1][0][0] + domino[1][1][0]
            line = '\t'.join([design_id] +
                             pair +
                             [_get_domino_id(seq, pair, domino_pairs)] +
                             [seq] +
                             [str(val) for sublist in list(domino[1])
                              for val in list(sublist)])

            overview_file.write(line + '\n')

        overview_file.write('\n')

    overview_file.close()

    return domino_pairs


def _write_dominoes_file(domino_pairs, filename):
    '''Writes an order file.'''
    order_file = open(filename + '.xls', 'w+')
    order_file.write(
        '\t'.join(['Well position', 'Name', 'Sequence', 'Notes']) + '\n')

    name_well_pos = {}
    pos = 0

    for domino, pairs in domino_pairs.iteritems():
        well_pos = plate_utils.get_well_pos(pos)
        name = '_'.join(pairs[0])
        name_well_pos[name] = well_pos
        pairs_str = ', '.join(set([pair[0] + '_' + pair[1] for pair in pairs]))
        order_file.write('\t'.join([well_pos, name, domino,
                                    pairs_str]) + '\n')
        pos += 1

    order_file.close()
    return name_well_pos


def _write_domino_pools_file(designs, domino_pairs, name_well_pos,
                             filename, vol=1.0):
    '''Writes domino pools operation file.'''
    design_id_well_pos = {}
    domino_file = open(filename + '_domino_pools.xls', 'w+')
    domino_file.write(
        '\t'.join(['DestinationPlateBarcode',
                   'DestinationPlateWell',
                   'SourcePlateBarcode',
                   'SourcePlateWell',
                   'ComponentName',
                   'Volume']) + '\n')

    pos = 0

    for design_id, design in designs.iteritems():
        for domino in design['dominoes']:
            well_pos = plate_utils.get_well_pos(pos)
            seq = domino[1][0][0] + domino[1][1][0]
            domino_id = _get_domino_id(seq, list(domino[0]), domino_pairs)
            design_id_well_pos[design_id] = well_pos

            line = '\t'.join([filename.split('/')[-1] + '_domino_pools',
                              well_pos,
                              filename.split('/')[-1],
                              name_well_pos[domino_id],
                              domino_id,
                              str(vol)])

            domino_file.write(line + '\n')

        pos += 1

    domino_file.close()
    return design_id_well_pos


def _write_parts_file(designs, filename):
    '''Writes a parts list.'''
    parts_well_pos = {}
    parts_file = open(filename + '_parts.xls', 'w+')
    parts_file.write(
        '\t'.join(['Well', 'Part / plasmid id']) + '\n')

    pos = 0
    for design in designs.values():
        for domino in design['dominoes']:
            for part in domino[0]:
                if part not in parts_well_pos:
                    well_pos = plate_utils.get_well_pos(pos)
                    parts_file.write('\t'.join([well_pos, part]) + '\n')
                    parts_well_pos[part] = well_pos
                    pos += 1

    parts_file.close()
    return parts_well_pos


def _write_plasmids_file(designs, design_id_well_pos, part_well_pos,
                         filename, vol=1.0):
    '''Writes plasmids file.'''
    plasmids_file = open(filename + '_plasmids.xls', 'w+')
    plasmids_file.write(
        '\t'.join(['DestinationPlateBarcode',
                   'DestinationPlateWell',
                   'SourcePlateBarcode',
                   'SourcePlateWell',
                   'ComponentName',
                   'Volume']) + '\n')

    pos = 0

    for design_id, design in designs.iteritems():
        well_pos = plate_utils.get_well_pos(pos)
        parts = []
        for domino in design['dominoes']:
            for part in list(domino[0]):
                if part not in parts:
                    line = '\t'.join([filename.split('/')[-1] + '_plasmids',
                                      well_pos,
                                      filename.split('/')[-1] + '_parts',
                                      part_well_pos[part],
                                      part,
                                      str(vol)])
                    plasmids_file.write(line + '\n')
                    parts.append(part)

        line = '\t'.join([filename.split('/')[-1] + '_plasmids',
                          well_pos,
                          filename.split('/')[-1] +
                          '_domino_pools',
                          design_id_well_pos[design_id],
                          design_id + '_domino_pool',
                          str(vol)])
        plasmids_file.write(line + '\n')

        pos += 1

    plasmids_file.close()


def _get_domino_id(seq, pair, domino_pairs):
    '''Gets a temporary domino id.'''
    if seq not in domino_pairs:
        domino_pairs[seq] = [pair]
    else:
        domino_pairs[seq].append(pair)

    return '_'.join(domino_pairs[seq][0])
