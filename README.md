# BIIN-Project

Disease Dynamics speeds up the research process of drug repurposing by automating literature searching and text mining in the pubmed database. This is powered by natural language processing techniques which is very efficient in searching  through large bodies of text. Automating text mining assists researchers in searching, filtering, and understanding hundreds of articles to identify patterns and develop connections from the text at a much faster pace. The user will be able to input biological mechanisms of the disease they wish to explore and they will get a list of relevant articles from pubmed and identified patterns of common medications in these groups. The intention of this tool is to assist researchers who wish to start researching about potential medications that can be repurposed.



My capstone project is made up of three files (1 python file and 2 html files) that can be downloaded to run the flask application over the local server at http://127.0.0.1:5000. 

It is very easy to run the application
  1. Download the three files included in this repository (BIIN_FINALPROJECT.py, index.html, search_results.html)
  2. Open the Command Prompt
  3. Set the directory to where the three files are now saved on your pc
  4. Type the command "python BIIN_FINAL_PROJECT.py" in the Command Prompt and click the enter button
  5. Wait a few seconds until you get the following message and then open the url on your browser

      * Serving Flask app 'BIIN_FINAL_PROJECT'
      * Debug mode: on
      WARNING: This is a development server. Do not use it in a production deployment. Use a production WSGI server instead.
      * Running on http://127.0.0.1:5000
      Press CTRL+C to quit
      * Restarting with stat
      * Debugger is active!
      * Debugger PIN: 135-337-769
    
You should now be on the homepage of the flask app. You can enter a couple keywords about the techniques used to find drugs, mechanisms of a disease, genes affected by a disease, and so on to explore relevant articles that will help start the accelerated research process of drug discovery. When you press the "Search" button you might have to wait a minute or two before getting any results. 

The results page will display a few things. The user inputted keywords will be displayed in the header, a list of medications found throughout the 100 most relevant articles, and the titles and abstract as well their relavancy scores of the 100 most relevant articles.

My report is in the file named BIIN Final Project Report.pdf and my presentation is in the file BIIN Final Project.pptx. I also included a link to my video demonstration of the tool working. All you have to do is make sure you sign in with a Ramapo email to view the link in Google Drive.

Hope this helps jumpstart your research towards discovering new uses for existing medications!
