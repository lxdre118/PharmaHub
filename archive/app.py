from flask import Flask, request, redirect, url_for, send_from_directory, render_template, jsonify
import os
from werkzeug.utils import secure_filename
from functions import allowed_file, process_structure_file, handle_file_upload, handle_file_download

UPLOAD_FOLDER = 'static/uploads'
DOWNLOAD_FOLDER = 'static/downloads'
ALLOWED_EXTENSIONS = {'sdf', 'smi', 'txt'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        return handle_file_upload(app, request)
    return render_template('upload.html')

@app.route('/waiting')
def waiting():
    return render_template('waiting.html')

@app.route('/result')
def result():
    return render_template('result.html')

@app.route('/downloads/<filename>')
def download_file(filename):
    return handle_file_download(app, filename)

if __name__ == '__main__':
    app.run(debug=True)
