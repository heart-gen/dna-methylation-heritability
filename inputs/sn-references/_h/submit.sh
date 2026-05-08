#!/bin/bash

## https://doi.org/10.1126/science.add7046
wget https://datasets.cellxgene.cziscience.com/0bb46187-dcf9-4b59-a825-33b6de1de21f.rds
mv -v 0bb46187-dcf9-4b59-a825-33b6de1de21f.rds caudate_body.rds

## https://doi.org/10.1038/s41467-024-45165-7
wget https://datasets.cellxgene.cziscience.com/cc8fa27f-1c5a-4ba8-bcf9-17b443fc6aca.rds
mv -v cc8fa27f-1c5a-4ba8-bcf9-17b443fc6aca.rds caudate_putamen.rds

## https://doi.org/10.1126/science.add7046
wget https://datasets.cellxgene.cziscience.com/9c504f0c-8332-413d-87ed-0d0f6ae5e290.rds
mv -v 9c504f0c-8332-413d-87ed-0d0f6ae5e290.rds hippocampal.rds

## https://doi.org/10.1126/science.add7046
wget https://datasets.cellxgene.cziscience.com/44fb7ba4-84e7-4073-89d0-dc6c2ce65782.rds 
mv -v 44fb7ba4-84e7-4073-89d0-dc6c2ce65782.rds dentate_gyrus.rds

## https://doi.org/10.1126/science.adf6812
wget https://datasets.cellxgene.cziscience.com/a2e7b6e9-79d5-446d-8fb1-19128c97a48f.rds
mv -v a2e7b6e9-79d5-446d-8fb1-19128c97a48f.rds dlpfc.rds
