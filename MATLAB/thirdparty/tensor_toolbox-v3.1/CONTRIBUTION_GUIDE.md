# Tensor Toolbox Contribution Guide

## Checklist

- [ ] **Issue** Before the merge request, submit an issue for the change, providing as much detailed information as possible. For bug reports, please provide enough information to reproduce the problem. 

- [ ] **Fork** Create a branch or fork of the code and make your changes.

- [ ] **Help Comments** Create or update comments for the m-files, following the style of the existing files. Be sure to explain all code options.

- [ ] **HTML Documentation** For any major new functionality, please follow the following steps.
  - [ ] Add HTML documentation in the `doc\html` directory with the name `XXX_doc.html`
  - [ ] Use the MATLAB `publish` command to create a new file in `doc\html` 
  - [ ] Add a pointer to this documentation file in `doc\html\helptoc.xml`
  - [ ] Add pointers in any related higher-level files, e.g., a new method for CP should be referenced in the `cp.html` file
  - [ ] Add link to HTML documentation from help comments in function
  - [ ] Update search database by running: builddocsearchdb('[full path to tensor_toolbox/doc/html directory]')
  
- [ ] **Tests** Create or update tests in the `tests` directory, especially for bug fixes or strongly encouraged for new code.

- [ ] **Contents** If new functions were added to a class, go to the `maintenance` directory and run `update_classlist('Class',XXX)` to add the new functions to the class XXX help information. If new functions were added at 
top level, go to `maintenance` and run `update_topcontents` to update the Contents.m file at the top level.

- [ ] **Release Notes** 
Update `RELEASE_NOTES.txt` with any significant bug fixes or additions.

- [ ] **Contributors List**
Update `CONTRIBUTORS.md` with your name and a brief description of the contributions.

- [ ] **Pass All Tests**
Confirm that all tests (including existing tests) pass in `tests` directory.

- [ ] **Merge Request** At any point, create a work-in-progress merge request, referencing the issue number and with this checklist and WIP in the header.


