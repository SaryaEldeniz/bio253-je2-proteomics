# Web-resource client offline checks

17 checks passed.

- PASS: measurement tables cannot be sent as identifiers
- PASS: empty and NA identifiers are rejected
- PASS: organism is explicit
- PASS: unresolved inputs are retained
- PASS: multiple returned candidates are labelled ambiguous
- PASS: wrong-species mappings require review
- PASS: duplicate input rows remain traceable
- PASS: service best match is not called exact
- PASS: mixed STRING organism namespaces are rejected
- PASS: UniProt exact accession and organism are retained
- PASS: redirected accession is distinguished
- PASS: unresolved UniProt lookup remains explicit
- PASS: KEGG requires explicit entries rather than a whole organism/database
- PASS: KEGG accepts organism genes and pathway entry identifiers
- PASS: KEGG headerless links retain the first data row and both columns
- PASS: KEGG no-link response returns a typed empty table
- PASS: Malformed KEGG responses are not silently parsed

No live API requests were made by these checks.
