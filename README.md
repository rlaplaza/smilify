# smilify
A wrapper for the xyz2mol flavor we use in cell2mol to turn an xyz into an rdkit mol with bonds. Putting this here because I somehow always end up having to write this again.
Note that smilify is significantly slower than the rdkit rdDetermineBonds implementatio and also uses the xyz2mol code as introduced by Jan Jensen as its core. However, we tried our best to extend it to other elements. On top of that, this one loops over covalent thresholds, which should make it a bit more robust in funky cases.
