# -*- coding: utf-8 -*-
# Generated by Django 1.11.4 on 2018-06-25 09:43
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='analysis_group',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=100, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='analysis_type',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('orderedFootprintList', models.CharField(max_length=200, null=True)),
                ('MajRefSeqSet', models.CharField(max_length=40, null=True)),
                ('listOfCladeCollectionsFoundInInitially', models.CharField(max_length=5000, null=True)),
                ('listOfCladeCollections', models.CharField(max_length=10000, null=True)),
                ('footprintSeqAbundances', models.CharField(max_length=10000, null=True)),
                ('footprintSeqRatios', models.CharField(max_length=100000, null=True)),
                ('clade', models.CharField(max_length=1)),
                ('coDom', models.BooleanField(default=False)),
                ('name', models.CharField(max_length=1000, null=True)),
                ('maxMinRatios', models.CharField(max_length=10000, null=True)),
                ('species', models.CharField(max_length=200, null=True)),
                ('artefactIntras', models.CharField(default='', max_length=5000)),
                ('isLockedType', models.BooleanField(default=False)),
                ('analysisGroupOf', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to='dbApp.analysis_group')),
            ],
        ),
        migrations.CreateModel(
            name='clade_collection',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('clade', models.CharField(max_length=1)),
                ('footPrint', models.CharField(default=True, max_length=100000)),
            ],
        ),
        migrations.CreateModel(
            name='clade_collection_type',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('analysisTypeOf', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.analysis_type')),
                ('cladeCollectionFoundIn', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.clade_collection')),
            ],
        ),
        migrations.CreateModel(
            name='data_analysis',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('listOfDataSubmissions', models.CharField(max_length=100, null=True)),
                ('withinCladeCutOff', models.FloatField(default=0.04)),
                ('typeSupport', models.FloatField(default=0.01)),
                ('cladeCollectionPopulationComplete', models.BooleanField(default=False)),
                ('analysisTypesDefined', models.BooleanField(default=False)),
                ('initialTypeDiscoComplete', models.BooleanField(default=False)),
                ('analysisTypesAssigned', models.BooleanField(default=False)),
                ('analysisTypesCollapsed', models.BooleanField(default=False)),
                ('refSeqsNamed', models.BooleanField(default=False)),
                ('speciesAssociated', models.BooleanField(default=False)),
                ('name', models.CharField(max_length=100, null=True)),
                ('description', models.CharField(max_length=5000, null=True)),
                ('dbVersionAnalysis', models.BooleanField(default=False)),
                ('dbVersionToRunAgainstID', models.IntegerField(null=True)),
                ('dSIDToAnalyse', models.IntegerField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='data_set',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='something', max_length=60)),
                ('reference_fasta_database_used', models.CharField(default='None', max_length=60)),
                ('submittingUser', models.CharField(default='Bob', max_length=30)),
                ('workingDirectory', models.CharField(default='None', max_length=300)),
                ('dataProcessed', models.BooleanField(default=False)),
                ('initialDataProcessed', models.BooleanField(default=False)),
                ('currentlyBeingProcessed', models.BooleanField(default=False)),
                ('timeStamp', models.CharField(default='None', max_length=40)),
            ],
        ),
        migrations.CreateModel(
            name='data_set_sample',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='None', max_length=200)),
                ('initialTotSeqNum', models.IntegerField(default=0)),
                ('post_seq_qc_absolute_num_seqs', models.IntegerField(default=0)),
                ('initialUniqueSeqNum', models.IntegerField(default=0)),
                ('finalTotSeqNum', models.IntegerField(default=0)),
                ('finalUniqueSeqNum', models.IntegerField(default=0)),
                ('non_sym_absolute_num_seqs', models.IntegerField(default=0)),
                ('nonSymSeqsNum', models.IntegerField(default=0)),
                ('size_violation_absolute', models.IntegerField(default=0)),
                ('size_violation_unique', models.IntegerField(default=0)),
                ('post_med_absolute', models.IntegerField(default=0)),
                ('post_med_unique', models.IntegerField(default=0)),
                ('initialProcessingComplete', models.BooleanField(default=False)),
                ('finalProcessingComplete', models.BooleanField(default=False)),
                ('errorInProcessing', models.BooleanField(default=False)),
                ('errorReason', models.CharField(default='noError', max_length=100)),
                ('cladalSeqTotals', models.CharField(max_length=5000, null=True)),
                ('dataSubmissionFrom', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.data_set')),
            ],
        ),
        migrations.CreateModel(
            name='data_set_sample_sequence',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('abundance', models.IntegerField(default=0)),
                ('cladeCollectionTwoFoundIn', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.clade_collection')),
                ('data_set_sample_from', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.data_set_sample')),
            ],
        ),
        migrations.CreateModel(
            name='reference_sequence',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='noName', max_length=30)),
                ('hasName', models.BooleanField(default=False)),
                ('clade', models.CharField(max_length=30)),
                ('sequence', models.CharField(max_length=500)),
                ('accession', models.CharField(max_length=50, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='symportal_framework',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('latest_reference_fasta', models.CharField(default='symClade_2_2.fa', max_length=30)),
                ('next_reference_fasta_iteration', models.IntegerField(default=1)),
                ('required_sub_e_value_seq_support_samples', models.IntegerField(default=3)),
                ('required_sub_e_value_seq_support_blast_symbiodinium', models.IntegerField(null=2)),
            ],
        ),
        migrations.AddField(
            model_name='data_set_sample_sequence',
            name='referenceSequenceOf',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.reference_sequence'),
        ),
        migrations.AddField(
            model_name='clade_collection',
            name='dataSetSampleFrom',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.data_set_sample'),
        ),
        migrations.AddField(
            model_name='analysis_type',
            name='dataAnalysisFrom',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.data_analysis'),
        ),
        migrations.AddField(
            model_name='analysis_group',
            name='dataAnalysisFrom',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='dbApp.data_analysis'),
        ),
    ]
