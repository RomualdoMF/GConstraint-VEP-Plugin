=head1 LICENSE

Copyright [2026] Romualdo Morandi Filho

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 romualdofarm@gmail.com
    
=cut

=head1 NAME

 GConstraint

=head1 SYNOPSIS

 mv GConstraint.pm ~/.vep/Plugins

 Annotate only the default columns (syn.z_score, syn.oe_ci.lower, syn.oe_ci.upper, mis.z_score, mis.oe_ci.lower, mis.oe_ci.upper, lof.z_score, lof.oe_ci.lower, lof.oe_ci.upper, lof.pLI):
 ./vep -i variations.vcf --plugin GConstraint,file=/path/to/gnomad/gnomad_constraint_final.tsv.gz

 Annotate any specified columns (e.g. syn.z_score, mis.z_score, lof.pLI):
 ./vep -i variations.vcf --plugin GConstraint,file=/path/to/gnomad/gnomad_constraint_final.tsv.gz,cols=syn.z_score,mis.z_score,lof.pLI

=head1 DESCRIPTION

 An Ensembl VEP plugin that adds gene constraint annotations from gnomAD version 4.1.1 (https://gnomad.broadinstitute.org/downloads#v4).

 Please cite the Mutational Constraint publication alongside the VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/

 The constraint metrics file can be downloaded from
 GRCh38: https://gnomad.broadinstitute.org/downloads#v4-constraint (Constraint metrics TSV)

 This file can be tabix-processed by:
 zcat gnomad.v4.1.1.constraint_metrics.tsv.bgz | (head -n 1 && tail -n +2  | sort -t$'\t' -k 9,9 -k 10,10n ) > gnomad.v4.1.1.constraint_metrics_sorted.tsv
 sed '1s/.*/#&/' gnomad.v4.1.1.constraint_metrics_sorted.tsv > gnomad.v4.1.1.constraint_metrics_final.tsv
 bgzip gnomad.v4.1.1.constraint_metrics_final.tsv
 tabix -f -s 9 -b 10 -e 11 gnomad.v4.1.1.constraint_metrics_final.tsv.gz

=cut

package GConstraint;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  die "ERROR: tabix not found in PATH\n" unless `which tabix 2>&1` =~ /tabix$/;
  die "ERROR: GConstraint file not provided or not found!\n"
    unless defined($param_hash->{file}) && -e $param_hash->{file};
  $self->add_file($param_hash->{file});

  # Default columns
  my @default_cols = qw(
    syn.z_score syn.oe_ci.lower syn.oe_ci.upper
    mis.z_score mis.oe_ci.lower mis.oe_ci.upper
    lof.z_score lof.oe_ci.lower lof.oe_ci.upper
    lof.pLI
  );

  # If cols= was passed, use those columns
  my @cols = @default_cols;
  if (defined $param_hash->{cols}) {
    @cols = split(/,/, $param_hash->{cols});
  }
  $self->{cols} = \@cols;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;
  my %info;

  # Basic descriptions (from gnomAD README)
  my %desc = (
    "syn.z_score"       => "Synonymous z-score",
    "syn.oe_ci.lower"   => "Synonymous observed/expected lower CI",
    "syn.oe_ci.upper"   => "Synonymous observed/expected upper CI",
    "mis.z_score"       => "Missense z-score",
    "mis.oe_ci.lower"   => "Missense observed/expected lower CI",
    "mis.oe_ci.upper"   => "Missense observed/expected upper CI",
    "lof.z_score"       => "Loss-of-function z-score",
    "lof.oe_ci.lower"   => "Loss-of-function observed/expected lower CI",
    "lof.oe_ci.upper"   => "Loss-of-function observed/expected upper CI",
    "lof.pLI"           => "Probability of being loss-of-function intolerant"
  );

  foreach my $c (@{$self->{cols}}) {
    $info{$c} = $desc{$c} || "Constraint annotation $c";
  }

  return \%info;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;
  my $transcript = $tva->transcript;
  my $start = $vf->{start};
  my $end   = $vf->{end};
  ($start, $end) = ($end, $start) if $start > $end;

  my ($res) = grep {
    $_->{transcript_id} eq $transcript->stable_id;
  } @{$self->get_data($vf->{chr}, $start, $end)};

  return $res ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;
  my @values = split /\t/, $line;

  my $transcript_id = $values[2];
  my %result;

  foreach my $c (@{$self->{cols}}) {
    # get index by name
    my $idx = $self->get_index_for_column($c);
    $result{$c} = $values[$idx] if defined $idx;
  }

  return {
    transcript_id => $transcript_id,
    result        => \%result
  };
}

sub get_index_for_column {
  my ($self, $col) = @_;

  my %map = (
    "syn.z_score"       => 25,
    "syn.oe_ci.lower"   => 26,
    "syn.oe_ci.upper"   => 27,
    "mis.z_score"       => 46,
    "mis.oe_ci.lower"   => 47,
    "mis.oe_ci.upper"   => 48,
    "lof.z_score"       => 91,
    "lof.oe_ci.lower"   => 92,
    "lof.oe_ci.upper"   => 93,
    "lof.pLI"           => 110,
    "gene_flags"        => 17,
    "constraint_flags"  => 18
  );

  return $map{$col};
}

1;