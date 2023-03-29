import csv

with open("PGGB Authors - Sheet1.tsv", newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')

    affiliation_dict, affiliation_dict_r = {}, {}
    affiliation_index = 1
    contribution_dict, contribution_dict_r = {}, {}
    contribution_index = 1
    authors = []
    authors_by_last_name = []
    funding = []
    
    latex_code = ''
    for r in reader:
        last_name = r['LastName'].strip()
        middle = r['MiddleInitials'].strip()
        first_name = r['FirstName'].strip()
        email = r['Email'].strip()
        if middle:
            first_name += ' ' + middle
        affiliation = r['Affiliations (Department, Organization, Street, City, Postcode, State, Country)'].strip()
        contribution = r['Contribution (TO FILL BETTER)'].strip()
        rank = float(r["Rank"])
        funding_string = r["Funding"]
        if funding_string != "":
            funding.append(funding_string.strip())
        authors.append((rank, last_name, first_name, affiliation, contribution, email))

    authors.sort()
    for rank, last_name, first_name, affiliations_string, contributions_string, email in authors:
        # Work out affliations
        affiliations = [i.strip() for i in affiliations_string.split(";")]
        affiliation_indices = []
        for affiliation in affiliations:
            if affiliation not in affiliation_dict:
                affiliation_dict[affiliation] = affiliation_index
                affiliation_dict_r[affiliation_index] = affiliation
                affiliation_index += 1
            affiliation_indices.append(str(affiliation_dict[affiliation]))

        # Work out contributions
        contributions = [i.strip() for i in contributions_string.split(",")]
        for contribution in contributions:
            if len(contribution) > 0:
                if contribution not in contribution_dict:
                    contribution_dict[contribution] = []
                contribution_dict[contribution].append(first_name + " " + last_name)

        latex_code += f"\\author[{','.join(affiliation_indices)}]{{\\fnm{{{first_name}}} \\sur{{{last_name}}}}}\\email{{{email}}}\n"

    for i in range(1, affiliation_index):
        latex_code += f"\\affil[{i}]{{{affiliation_dict_r[i]}}}\n"

    print(latex_code)