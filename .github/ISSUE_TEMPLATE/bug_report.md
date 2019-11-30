---
name: Bug report
about: Something went wrong? Fill a bug report!
title: "[BUG] "
labels: ''
assignees: sergiocorreia

---

# Before submitting the bug report

1. Verify that you are using the latest versions of reghdfe and ftools (`which reghdfe`). Note that the latest version is usually on Github and not on SSC.
2. Verify that your Stata is updated (`update query`).

# Bug report

**Please complete the following information:**

- Stata version: [e.g. Stata 16.0 31dec2019] Find out with the `update` and `version` commands
- OS: [e.g. Windows 10, Linux 18.04]

**Behavior**

- Expected behavior: a clear and concise description of what you expected to happen.
- Actual behavior: what was the error. If possible, copy-paste output. For instance:

```stata
. reghdfe price weight length, absorb(foreign#)
nothing found where name expected
                 stata():  3598  Stata returned error
         fixed_effects():     -  function returned error
                 <istmt>:     -  function returned error
r(3598);
```

**Steps to reproduce the problem**

If possible, please provide a [MWE](https://en.wikipedia.org/wiki/Minimal_working_example) so I can reproduce the problem.
