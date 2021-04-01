# The Pyrfume Public Data Archive
## Wrangled, pre-processed, and curated olfactory psychophysics data.

![pyrfume-logo](https://avatars.githubusercontent.com/u/34174393)

Intended for use with the `pyrfume` Python library:
```console
pip install pyrfume
```

Examples:
```python
import pyrfume
```
followed by e.g.:

```python
# Load all the data from Bushdid et al, 2014 ("Humans Can Discriminate More than 1 Trillion Olfactory Stimuli")
molecules = pyrfume.load_data('bushdid_2014/molecules.csv', remote=True)
mixtures = pyrfume.load_data('bushdid_2014/mixtures.csv', remote=True)
behavior = pyrfume.load_data('bushdid_2014/behavior.csv', remote=True)
```
or e.g.:

```python
# Load all the data from Snitz et al, 2013 ("Predicting Odor Perceptual Similarity from Odor Structure")
molecules = pyrfume.load_data('snitz_2013/molecules-info.csv', remote=True)
behavior = pyrfume.load_data('snitz_2013/behavior-main.csv', remote=True)
```

Analogous examples can be found [here](code_examples.py) for other datasets.