# Contributing Guidelines  

Hi! It's great to have you here! In this project we will explore how ecological networks in bromeliads difer between locations. It's heavily inspired on [this paper](https://www.cambridge.org/core/journals/parasitology/article/abs/beta-and-phylogenetic-diversities-tell-complementary-stories-about-ecological-networks-biogeography/2FBB62E55AAD5EE67072C14D686BF72D#article), so maybe it's a good idea to start by taking a look on it.

We will use this repository to version control a manuscript, including code and figures. Please read this document carefully if you want to make a contribution.

## Code of Conduct

If you decide to participate in this project, we assume you comply with the Bromeliad Working Group [Constitution](https://www.zoology.ubc.ca/~srivast/bwg/constitution.html). By participating, you are expected to uphold this code. Please report unacceptable behavior to [graciellehigino@gmail.com](mailto:graciellehigino@gmail.com).  

## Before you start  

This is a project to investigate how ecological networks vary in space. We will use Julia to analyse data and create figures, and markdown to write our manuscript. You can read more on our questions and objectives in our [README](README.md) file, and find the taks we are working on in our [Roadmap](Roadmap.md).  
This is a completly open project, and every contribution made to it is also completly open. Keep that in mind, and let's move foward!  

### Project organization  

```shell
ms/  
    ⊢manuscript.md
    ⊢references.bib
    ⊢style.csl
    ⊢template.tex
code/
data/
figures/
tables/
.gitignore
LICENSE.md
Manifest.toml
Project.toml
README.md
CONTRIBUTING.md 
ROADMAP.md
```
 
## How to work on GitHub  

If you're new to GitHub, take some time to understand how a repository works. Don't worry if you find it difficult, we can always help!  
GitHub has its own [guides](https://guides.github.com/) which are really useful. Basically, you need to know how to [file an Issue](https://guides.github.com/features/issues/) and how to contribute directly to the core contents of the project by [submitting a pull request](https://guides.github.com/activities/forking/#making-a-pull-request). Please, be as clear as you can in both cases, describing what is the problem/suggestion/solution in details. Also, be sure that your suggestion/problem/solution was not addressed in another Issue before.

Each part of this project should have a branch (e.g. ms_intro, ms_methods, code_data, code_analysis1), and your contributions must be made in a branch with your initials (e.g. ms_intro_gh, code_data_gh). Please only make modifications on your branch unless clearly stated otherwise. When finished, make a pull request to merge your branch with the main branch of your section (e.g. ms_intro_gh -> ms_intro) and tag a reviewer to approve the merge.

**ALWAYS** make sure to pull changes before you start working, commit and stage frequently, and pull often - especially before breaks.

Thank you again for contributing! **HAVE FUN!**
