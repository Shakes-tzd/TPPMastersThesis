---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---


# Introduction

Mercury (Hg) is a severe health hazard in the environment and for people, especially fetuses and young children {cite}`gibb_mercury_2014`. The release of mercury into the atmosphere results from human activities and natural processes. Upon release, mercury moves between the air, soils, and waters until, eventually, it is removed from the system through burial in coastal and deep ocean sediments, lake sediments, and subsurface soils {cite}`esdaile_mercury_2018`. Moreover, mercury can travel great distances when emitted into the atmosphere, contaminating ecosystems, fish, birds, mammals, and human food chains {cite}`esdaile_mercury_2018`. Recent estimates, as seen in {numref}`gma2018_hg-emissions_by-industry`,  indicate that about 38\% of global Hg emissions come from ASGM, making it the largest source of Hg pollution to the atmosphere and hydrosphere around the world{cite}`united_nations_environment_programme_technical_2019`. However, the amount of Hg released by ASGM activities and the extent to which it is transported regionally and globally is highly uncertain. This thesis examines how top-down emissions estimation techniques can help reduce the uncertainty in ASGM Hg estimates. First, we will briefly describe the Hg problem in ASGM. As a next step, we review the state of atmospheric Hg monitoring and modeling in regions where ASGM activities are prevalent. In addition, we use observed atmospheric Hg concentrations and atmospheric Hg predictions produced by the GEOS-Chem Model to come up with the first-ever top-down estimates of ASGM Hg emissions for a region in Peru. Lastly, we synthesize the findings and make recommendations for policymakers to use monitoring and modeling to better understand ASGM Hg emissions. To comply with the Minamata Convention, we will use the term "Hg emissions" throughout this thesis to refer to Hg discharges into the atmosphere and "Hg releases" to refer to Hg discharges into the ground and water. "Hg discharge" will be used when referring to both emissions and releases.

```{figure} ../figures/07-24-22_gma2018_hg-emissions_by-industry.pdf
---
height: 200px
name: gma2018_hg-emissions_by-industry
---
Pie chart showing the 2018 global mercury assessment (GMA 2018) ASGM Hg emission estimates for different sectors. ASGM is estimated to emit the most Hg emissions (shown in orange) at 838 Mg, followed by industry sectors (shown in red) at 614 Mg, then fuel combustion (shown in blue) at 533 Mg, and finally, intentional use sectors excluding ASGM (show in green) at 239 Mg {cite}`united_nations_environment_programme_technical_2019`.
```

## Motivation

More than 100 million people depend on artisanal and small-scale gold mining (ASGM) for their livelihood globally, particularly in the over 81 countries, predominantly in the global south where ASGM exists{cite}`planetgold_planetgold_2021`. Additionally, ASGM is an essential source of income and an opportunity for rural development in countries where options and alternatives to ASGM for generating income to buy necessities of daily life are in short supply or nonexistent {cite}`planetgold_planetgold_2021`. It is estimated that around 10 to 20 million (ASGM) miners are employed in ASGM worldwide - about a third of them are women - and they provide 90\% of the global gold mining workforce and extract about 20\% of the world's gold annually {cite}`planetgold_planetgold_2021`. For example, in Peru, ASGM sustains the livelihoods of an estimated 1 million people, and between 300,000 and 500,000 miners were involved in Peru's ASGM sector as of 2014. Despite being a vital source of livelihood for the communities that practice ASGM, its activities often lead to several environmental, human, and social harms. In addition to Hg releases to the environment, ASGM externalities include deforestation, tropical diseases such as malaria, dangerous and unsafe working conditions, crime and exploitation of indigenous communities, diesel and gasoline spills, and human trafficking {cite}`usaid_usaid_2020`.

While most of the Hg pollution from ASGM is local, its ability to travel across borders and contaminate distant ecosystems bolsters the case for concerted global efforts to eliminate Hg pollution in all forms, including ASGM Hg pollution. The Minamata Convention (MC) is one of the unified global efforts to combat Hg pollution. It is a legally binding global treaty with 137 parties as of this writing, and its goal is to protect human health and the environment from the adverse effects of mercury. The MC's text comprises articles that address different sources of Hg pollution, including ASGM Hg emissions{cite}`unep_minamata_2013`. This thesis work is inspired by the need to improve our understanding of Hg emissions and Hg's regional and global transport to ensure that policies and actions taken to reduce ASGM Hg emissions are informed by the best available science and take advantage of all resources at our disposal to create actionable scientific knowledge.

## What are roles and directives?

Roles and directives are two of the most powerful tools in Jupyter Book. They
are kind of like functions, but written in a markup language. They both
serve a similar purpose, but **roles are written in one line**, whereas
**directives span many lines**. They both accept different kinds of inputs,
and what they do with those inputs depends on the specific role or directive
that is being called.

### Using a directive

At its simplest, you can insert a directive into your book's content like so:

````
```{mydirectivename}
My directive content
```
````

This will only work if a directive with name `mydirectivename` already exists
(which it doesn't). There are many pre-defined directives associated with
Jupyter Book. For example, to insert a note box into your content, you can
use the following directive:

````
```{note}
Here is a note
```
````

This results in:

```{note}
Here is a note
```

In your built book.

For more information on writing directives, see the
[MyST documentation](https://myst-parser.readthedocs.io/).

### Using a role

Roles are very similar to directives, but they are less-complex and written
entirely on one line. You can insert a role into your book's content with
this pattern:

```
Some content {rolename}`and here is my role's content!`
```

Again, roles will only work if `rolename` is a valid role's name. For example,
the `doc` role can be used to refer to another page in your book. You can
refer directly to another page by its relative path. For example, the
role syntax `` {doc}`index` `` will result in: {doc}`index`.

For more information on writing roles, see the
[MyST documentation](https://myst-parser.readthedocs.io/).

### Adding a citation

You can also cite references that are stored in a `bibtex` file. For example,
the following syntax: `` {cite}`holdgraf_evidence_2014` `` will render like
this: {cite}`holdgraf_evidence_2014`.

Multiple citations can be used like this:
 {cite}`holdgraf_rapid_2016, holdgraf_encoding_2017`

Moreover, you can insert a bibliography into your page with this syntax:
The `{bibliography}` directive must be used for all the `{cite}` roles to
render properly.
For example, if the references for your book are stored in `references.bib`,
then the bibliography is inserted with:

````
```{bibliography}
```
````

Resulting in a rendered bibliography that looks like:

<!-- ```{bibliography}

``` -->

### Executing code in your markdown files

If you'd like to include computational content inside these markdown files,
you can use MyST Markdown to define cells that will be executed when your
book is built. Jupyter Book uses _jupytext_ to do this.

First, add Jupytext metadata to the file. For example, to add Jupytext metadata
to this markdown page, run this command:

```
jupyter-book myst init markdown.md
```

Once a markdown file has Jupytext metadata in it, you can add the following
directive to run the code at build time:

````
```{code-cell}
print("Here is some code to execute")
```
````

When your book is built, the contents of any `{code-cell}` blocks will be
executed with your default Jupyter kernel, and their outputs will be displayed
in-line with the rest of your content.

For more information about executing computational content with Jupyter Book,
see [The MyST-NB documentation](https://myst-nb.readthedocs.io/).
