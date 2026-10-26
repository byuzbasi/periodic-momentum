"""Read raw CSV and XLSX feature cells without spreadsheet dependencies."""
from pathlib import Path
import csv
import json
import posixpath
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parent
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def column_number(cell_reference):
    value = 0
    for letter in cell_reference:
        if not letter.isalpha():
            break
        value = 26 * value + ord(letter.upper()) - ord('A') + 1
    return value - 1


def sheets(path):
    with zipfile.ZipFile(path) as z:
        rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        targets = {x.attrib['Id']: x.attrib['Target'] for x in rels}
        result = []
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).find('s:sheets', NS):
            target = targets[sheet.attrib[REL]]
            target = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)
            assert target.startswith('xl/') and '..' not in target.split('/')
            result.append((sheet.attrib['name'], target))
        return result


def xlsx_rows(path, sheet_name):
    target = dict(sheets(path))[sheet_name]
    with zipfile.ZipFile(path) as z:
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            strings = [''.join(t.text or '' for t in x.iterfind('.//s:t', NS))
                       for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        data = ET.fromstring(z.read(target)).find('s:sheetData', NS)
        result = []
        for row in data:
            cells = {}
            for c in row.findall('s:c', NS):
                typ = c.attrib.get('t', 'n')
                value = c.find('s:v', NS)
                value = value.text if value is not None else None
                if typ == 's' and value is not None:
                    value = strings[int(value)]
                elif typ == 'inlineStr':
                    value = ''.join(x.text or '' for x in c.iterfind('.//s:t', NS))
                cells[column_number(c.attrib['r'])] = (value, typ, c.find('s:f', NS) is not None)
            result.append((int(row.attrib['r']), cells))
        return result
