#! /bin/bash

pbpaste | sed "s/$1/$2/g" | pbcopy
