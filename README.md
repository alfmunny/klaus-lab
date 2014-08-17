## Foreword

This page contains the full steps and descriptions to generate models and results.

Results are represented in diagram above. Now there are still results for 12 dimensions left.
They would be filled in couple of days.

Only to see the results, please review them on [Github](https://github.com/alfmunny/klaus-lab).

For each model is a evaluation file `SUM` generated.
For example, the `SUM` for the `bakis` of `10dim-Slant`
can be found [here on Github](https://github.com/alfmunny/klaus-lab/blob/master/hmm/10dim-Slant/bakis/SUM).
The content of `SUM` will be described here in part **Evaluation of models**.

## Project Folders

Here are relevant folders:

`def`: definitions for lexicon

`ref`: character reference

`lm`: language model

`lin`: current line files

`lin_noSlant`: line image files without slant

`lin_Slant`: line image files with slant

`ufv`: current ufv (binary date of features) files

`ufv_10dim_noSlant`: features for 10 dimensions without slant

`ufv_12dim_noSlant`: features for 12 dims without slant

`ufv_10dim_Slant`: features for 10 dims with slant

`ufv_12dim_Slant`: features for 12 dims with slant

`hmm`: hidden markov model



## Configuration before run

It's also possible to generate the result by self.

Configuration to do:

* Download the .zip file and override the original folder with new folder in `lab/klaus/`.

* `doc/` folder is not including. Make sure you have it locally.

* override the `esmeralda/src/pen/lib/segmentation/fextract.c` with the `fextract.c` which provided in this project root
  and then compile

      $ cd esmeralda/src/pen
      $ make
      $ make install

## Generate lines

Direct to corresponding folder `lin_noSlant` or `lin_Slant` and run:

```
$ ./MakeLines
```

## Generate features

Change the current `lin` to the lines files which you want to generate from.

For example  when generate features with lines with slant, override the current `lin` with `lin_Slant`:

```
$ rm -r lin
$ cp -r lin_Slant lin

```

Then direct to the ufv folder with the desired configuration showed in the folder name.

Example:

```
$ cd ufv_10dim_Slant
$ ./merkmale_hmm.sh
```
or

```
$ cd ufv_12dim_Slant
$ ./merkmale_hmm.sh
```

When generate features without slant, override the current `lin` with `lin_noSlant` just like before:

```
$ rm -r lin
$ cp -r lin_noSlant lin

```
and run the `merkmale_hmm.sh` in `ufv_10dim_noSlant/` and`ufv_12dim_noSlant`.

## Generate models

In `hmm/` are these folders for storing the generated models.

* `10dim-Slant`: 10 features with slant
* `12dim-Slant`: 12 features with slant
* `10dim-noSlant`: 10 features without slant
* `12dim-noSlant`: 12 features without slant

Folder tree of `10dim-Slant`:

    10dim-Slant/
      | linear/  :initial training
      | nobakis/  :final linear(no bakis) models
      | bakis/    :final bakis models
      | DoAll    :shell script for generating models

**`DoAll` is a shell script, which run all commands from the very beginning
to the end.** It means it will override all current models, generating it
again and also do the evaluation for models in `bakis` and `nobakis`.

**Attention: please choose the corresponding ufv folder before generating the model.**

Example: `10dim-Slant`

Generate models for `10dim-Slant`, firstly make sure the current ufv folder is including the corresponding
features. You can override the current ufv folder with the features you have generated before.

```
$ rm -r ufv
$ cp -r ufv_10dim_Slant ufv
```

Then run `DoAll` in the folder with the corresponding name

```
$ cd hmm/10dim_Slant
$ ./DoAll
```

## Evaluation of models

### Introduction of Evaluation

After `DoAll` the evaluation summary `SUM` will be created in `bakis` and `nobakis`.

The evaluation format:

    WA/SR/WC/S/D/I = 60.85+/-1.2 0.00 65.65 29.11 4.61 4.79

Criteria:

`WA`: word accuracy

`SR`: sentence recognition

`WC`: words correct

`S`: substitutions

`D`: deletions

`I`: insertions

`WC` + `S` + `D` = 100%

`WA` = 100% - `S` - `D` - `I`

### Run evaluation separately

To only run evaluation saperately, direct to the model folder `bakis` or `nobakis`.

Here is a example for `10dim-Slant`

```
$ cd 10dim-Slant
$ cd bakis
$ ../../MkEvals
$ ../../MkSummary
```
The results will be printed directly in terminal.

Or stored the results in file `SUM`.

```
$ ../../MkSummary>SUM
```
