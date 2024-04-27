from flask import Flask, render_template, request
from Bio import Entrez
import xml.etree.ElementTree as ET
from nltk.tokenize import wordpunct_tokenize
from drug_named_entity_recognition import find_drugs
import re
from collections import Counter

app = Flask(__name__)

Entrez.email = 'eidil@ramapo.edu'  

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    search_terms = request.form.get('search_terms', '').split(',')
    search_query = ' AND '.join(search_terms)

    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=100)
    record = Entrez.read(handle)

    id_list = record['IdList']
    records, common_medication_names = get_pubmed_records(id_list, search_terms)  
    records.sort(key=lambda x: x['relevance_score'], reverse=True)

    return render_template('search_results.html', search_terms=search_terms, records=records, common_medication_names=common_medication_names)

def get_pubmed_records(id_list, search_terms):
    results = []
    common_medication_names = []  
    for pubmed_id in id_list:
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        xml_data = handle.read()
        root = ET.fromstring(xml_data)

        for article in root.findall('.//PubmedArticle/MedlineCitation/Article'):
            title = article.find('.//ArticleTitle').text if article.find('.//ArticleTitle') is not None else 'N/A'
            abstract = article.find('.//AbstractText').text if article.find('.//AbstractText') is not None else 'N/A'
           
            if title is None:
                title = 'N/A'
            if abstract is None:
                abstract = 'N/A'

            title_tokens = wordpunct_tokenize(title)
            abstract_tokens = wordpunct_tokenize(abstract)
            title_drugs = find_drugs(title_tokens, is_ignore_case=True)
            abstract_drugs = find_drugs(abstract_tokens, is_ignore_case=True)
            all_drugs = title_drugs + abstract_drugs

            pattern = r"'name':\s*'([^']+)'"
            medicine_names = []
            for drug in all_drugs:
                match = re.search(pattern, str(drug))
                if match:
                    medicine_name = match.group(1)
                    medicine_names.append(medicine_name)
                    
            medicine_name_counts = Counter(medicine_names)
            common_medication_names = list(medicine_name_counts.items())

            relevance_score = calculate_relevance_score(search_terms, title, abstract)

            results.append({
                'title': title,
                'abstract': abstract,
                'relevance_score': relevance_score  
            })

    return results, common_medication_names

def calculate_relevance_score(search_terms, title, abstract):
    score = 0

    for term in search_terms:
        score += title.lower().count(term.lower())
        score += abstract.lower().count(term.lower())
    return score

if __name__ == '__main__':
    app.run(debug=True)


