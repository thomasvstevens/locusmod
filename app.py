"""Flask app serving submission form and download pages"""

from flask import Flask, g, render_template, request, send_file
import batch
import locus
import os
import shutil

app = Flask(__name__)
files = ['primers.csv', 'idt_import.csv', 'plasmids.zip', 'loci.zip']
marker_keys = [k for k in batch.MAP_PLASMID_PATH.keys() if k.startswith('Lox')]
index_kw = {'mod_keys':batch.MAP_MOD_TYPE.keys(),
            'assembly_keys':batch.MAP_ASSEMBLY_METHOD.keys(),
            'marker_keys':marker_keys}

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', **index_kw)

@app.route('/download', methods=['GET', 'POST'])
def download():
    if request.method == 'POST':
        batch_kw = {k:v for k,v in request.form.items()
                if not k.startswith('start')}
        b = batch.Batch(**batch_kw)
        shutil.rmtree(batch.paths.OUTPUT_PREFIX)
        os.mkdir(batch.paths.OUTPUT_PREFIX)
        b.write_primers_csv()
        b.assign_numbers(request.form['start_oligos'])
        b.write_idt_csv()
        b.write_plasmids_zip(request.form['start_plasmids'])
        b.write_loci_zip()
        rows = b.list_operations()
        return render_template('download.html', files=files, rows=rows)
    return render_template('index.html', **index_kw)

def test():
    with app.test_request_context('/', method='GET'):
        assert request.path == '/'
        assert request.method == 'GET'
    with app.test_request_context('/', method='POST'):
        assert request.path == '/download'
        assert request.method == 'POST'

if __name__=='__main__':
    app.run()
    test()
