# The Goodscents Company Dataset

This is a CSV dump of a web scraped HTML data of [The Goodscents Company website by William Luebke](http://www.thegoodscentscompany.com/). The web scraping and the CSV transformation is maintained by [Contrebande Labs](https://github.com/contrebande-labs/) in [their fork](https://github.com/contrebande-labs/pyrfume-data) of the [Pyrfume Data repository](https://github.com/pyrfume/pyrfume-data). This dataset is an extraction of a web scraping from January 17 2022. The data has not been modified, transformed or filtered except when specifically indicated.

## Usage

### `data_rw_odor.csv`

You can get a list of raw (rw) aromachemical IDs and their corresponding olfactive descriptors from the `data_rw_odor.csv` file. The oranoleptic tags are in the `Tags` column and the Goodscents raw aromachemical profile ID are in the `TGSC ID` column. For some profiles the organoleptic data was provided as a full english sentence, in these cases, some slight parsing was applied to obtain an array of descriptor strings. When the `Original Tags` column is set to `false`, the organoleptic tag array in the `Tags` column has been derived from the plain text value in the `Description` column, otherwise the tag array is left intact from the HTML content.

### `data_rw_contents.csv` (optional)

Optionally, you can link more raw aromachemical profile IDs to the previously obtained organoleptic data by looking at the `data_rw_contents.csv`. From that file, in the `TGSC ID` column you will find the "mixture/solution" raw (rw) aromachemical profile IDs and in the `Content TGSC ID` the profile ID for products/aromachemicals that might be included in the "mixture/solution". These "contained" products are raw aromachemicals profile when the value in the `Content Data Category` is set to "rw". Beware that some "mixture/solution" contents are actually solvents or impurities and in any case may not share the same organoleptic profile as their parent "mixture/solution".

### `data_rw_opl.csv`

Then you need the molecular structure(s) profile (OPL) ID(s) associated with each raw aromachemical profile ID obtained previously by looking at the `data_rw_opl.csv`. In this file, you will find the molecular structure profile OPL ID in the `TGSC OPL ID` column and the raw aromachemical ID it is associated with in the `TGSC ID` column.

### `opl.csv`

Finally, the molecular structure data (SMILES, InChI, etc.) can be found in the `opl.csv` file, indexed by an OPL ID (`TGSC OPL ID`) previously obtained from the `data_rw_opl.csv` file.

## Known Issues
At the moment, only the raw (rw) aromachemical molecular structure (SMILES, InChi, etc.) to olfactive activity (descriptor labels) data is provided. We will soon integrate the following datapoints :

* Flavor and fragrance demo formulas (i.e. olfactive activity of mixtures of raw aromachemicals);
* Natural products (essential oils, absolutes, concretes, tinctures, extractions, etc.) and their typical GC/MS profiles (natural occurence/provenance/source of raw aromachemicals);
* Accessory data points from the raw (rw) aromachemical profiles (PubChem, IFRA, etc.).