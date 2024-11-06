import camelot

def extract_tables_with_camelot(pdf_path):
    tables = camelot.read_pdf(pdf_path, pages="all")
    return tables

tables = extract_tables_with_camelot("./crust/crust_Hebeler_2013_ApJ_773_11.pdf")

print(tables)