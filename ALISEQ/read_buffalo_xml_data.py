# -*- coding: utf-8 -*-

"""

.. Created on Thu Aug  9 16:45:06 2018

   @author: cow082

   python read_buffalo_xml_data.py

"""

def processXML():
    '''used to construct 'resources/buffalo_sample_metadata.csv'
    
    # it uses the following xml files 
    "ERS2495773-ERS2495833.xml", "ERS2495835-ERS2495900.xml","ERS2495902-ERS2495992.xml"
    
    # which are obtained from here (there is an xml download link on the page): 
    https://www.ebi.ac.uk/ena/data/view/ERS2495773-ERS2495833
    https://www.ebi.ac.uk/ena/data/view/ERS2495835-ERS2495900
    https://www.ebi.ac.uk/ena/data/view/ERS2495902-ERS2495992
    '''
    
    PATH = "C:/Users/cow082/Desktop/roslin secondment documentation"
    
    infiles = [
        "ERS2495902-ERS2495992.xml",
        "ERS2495835-ERS2495900.xml",
        "ERS2495773-ERS2495833.xml",
        ]
    
    outlines = []
    for f in infiles:
        contents = open("%s/%s" % (PATH, f), 'r').read()
        chunks = contents.split('<SAMPLE ')[1:]
        for chunk in chunks:
            lines = chunk.split('\n')
            for line in lines:
                if '<PRIMARY_ID>' in line:
                    name = line.split('<PRIMARY_ID>')[1].split('</PRIMARY_ID>')[0]
                elif '<TITLE>' in line:
                    data = line.split('<TITLE>')[1].split('</TITLE>')[0]
                    components = data.split(',')
                    data = [item.strip() for item in ','.join(components[:-3]).split(' from an ') + components[-3:]]
                    data = ','.join(data)
            outlines.append(name+','+data)
            
    
    with open('%s/buffalo_sample_metadata.csv' % PATH, 'w') as F:
        F.write('\n'.join(outlines))


def main():
    processXML()
    
    
if __name__ == '__main__':
    main()

