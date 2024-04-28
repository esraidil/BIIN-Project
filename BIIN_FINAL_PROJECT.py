from flask import Flask, render_template, request
from Bio import Entrez
import xml.etree.ElementTree as ET
from nltk.tokenize import wordpunct_tokenize
from drug_named_entity_recognition import find_drugs
import re

app = Flask(__name__)

Entrez.email = 'eidil@ramapo.edu'  

#This section runs the index.html file which holds the home page of the web app. It is also where the user input is entered.
@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

#This section runs the search_results.html file which holds the results page of the web app. It is also where the results returned from the user input is displayed.
@app.route('/search', methods=['POST'])
def search():

    #User enters keywords that are then split by the delimeter comma and stored separately as search_terms then joined to form search_query.
    search_terms = request.form.get('search_terms', '').split(',') 
    search_query = ' AND '.join(search_terms) 

    #Entrez is used to search the PubMed database for the search_query words and returnes the 100 most relevant articles.
    handle = Entrez.esearch(db="pubmed", term=search_query, retmax=100)
    record = Entrez.read(handle)

    #This section runs the function pubmed_records for the two main types of results, records and common_med_names.
    ID_List = record['IdList']
    records, common_med_names = pubmed_records(ID_List, search_terms)  
    records.sort(key=lambda x: x['relevance_score'], reverse=True) #Sorts the column of the table named relevance_score in descending order.
    common_med_names = sorted(common_med_names)

    #Finally, the search_results.html file is rendered and is given inofrmation from the search_terms, records, and common_med_names.
    return render_template('search_results.html', search_terms=search_terms, records=records, common_med_names=common_med_names)

#The relevance_score is calculated for each row in the table by counting the number of occurrances of all search terms in the title and abstract.
def calculate_relevance_score(search_terms, title, abstract):
    score = 0
    for term in search_terms:
        score += title.lower().count(term.lower())
        score += abstract.lower().count(term.lower())
    return score

#The function used to pull pubmed articles from the database and further pull information such as title and abstract from these records.
def pubmed_records(ID_List, search_terms):
    table_results = [] #This list holds the title, abstract, and rlevance_score information found from pubmed articles. 
    common_med_names = [] #This list holds the common medication names that show up more than once in the title and abstract returned by table_results.

    #The for loop goes through each pubmed article retrieved from pubmed.
    for pubmed_ID in ID_List:
        handle = Entrez.efetch(db="pubmed", id=pubmed_ID, retmode="xml")
        xml_data = handle.read()
        root = ET.fromstring(xml_data)

        #The for loop connects to the three web pages and retrieves the data as article, title, and abstract.
        for article in root.findall('.//PubmedArticle/MedlineCitation/Article'):
            title = article.find('.//ArticleTitle').text if article.find('.//ArticleTitle') is not None else 'N/A'
            abstract = article.find('.//AbstractText').text if article.find('.//AbstractText') is not None else 'N/A'

            #Makes sure that null data is replaced with 'N/A' to avoid errors.
            if title is None:
                title = 'N/A'
            if abstract is None:
                abstract = 'N/A'

            #The wordpunct_tokenize function from the NLTK tokenizer package separates words in a sentence so that each word can be analyzed separately. 
            title_tokens = wordpunct_tokenize(title)
            abstract_tokens = wordpunct_tokenize(abstract)

            #The find_drugs function from the drug_named_entity_recognition package compares the list of words from the tokens output against the list of available drugs in their database and flags matches.
            title_drugs = find_drugs(title_tokens, is_ignore_case=True) 
            abstract_drugs = find_drugs(abstract_tokens, is_ignore_case=True)
            all_drugs = title_drugs + abstract_drugs

            pattern = r"'name':\s*'([^']+)'" #The raw data obtained from the find_drugs function is in the form of the drug name and all of its synonyms so the pattern specifies that the program is only interested in the drug name.
            
            #This for loop uses the re module to form a match between the pattern and string of drugs found, if there is a match, it will place it in the medicine_names list.
            for drug in all_drugs:
                match = re.search(pattern, str(drug))
                if match:
                    medicine_name = match.group(1)
                    common_med_names.append(medicine_name)
            common_med_names = list(set(common_med_names))

            #The relevance_score function is used to assign a score for the combination of the title and abstract.
            relevance_score = calculate_relevance_score(search_terms, title, abstract)

            #The table is created for the results from articles. 
            table_results.append({
                'title': title,
                'abstract': abstract,
                'relevance_score': relevance_score  
            })

    return table_results, common_med_names

if __name__ == '__main__':
    app.run(debug=True)


