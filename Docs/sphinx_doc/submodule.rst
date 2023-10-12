.. _Submodule:

Git Submodules
--------------

When using the submodule to build, it is ideal to properly update and match what is in the repository. Depending on Git version, different commands and options to ensure these match. An example workflow is to run ``git pull`` to get the latest commit on your current branch, and then run ``git submodule update`` to explicitly update the submodule. This should work for all versions of ``git`` which support submodules.

The following example demonstrates a shorter form that combines both commands and requires Git 2.14 or newer:

    .. code:: shell
              # Replaces your git pull to use both the updated code and the updated submodule
              git pull --recurse-submodules
The following example demonstrates setting defaults in the config file for the repository and requires Git 2.15 or newer:

    .. code:: shell
              # Set this once
              git config submodule.recurse true
              # When configured as true, this will use both the updated code and the updated submodule
              git pull
              # This option will override any configuration and not update the submodule
              git pull --no-recurse-submodules
These example also apply to ``git checkout``. For more details see _`Git Tools Submodules`: https://git-scm.com/book/en/v2/Git-Tools-Submodules
