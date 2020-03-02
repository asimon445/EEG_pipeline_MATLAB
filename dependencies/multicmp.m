%  Function File: multicmp
%
%          USAGE: [padj,S] = multicmp (y,g,type,pairs,option)
%
%             or  [padj] = multicmp (p,option)
%
%
%  A set of multiple comparison tools that apply procedures to control
%  either the family-wise error rate (FWER) or the false discovery rate
%  (FDR) when performing post-tests. The FWER is the probability of
%  making at least one false discovery (type 1 error) among the entire
%  set of hypotheses when multiple testing. In contrast, the FDR is the
%  expected proportion of false discoveries as a fraction the number of
%  null hypothesis rejected among a set of hypotheses. For example, for
%  an experimental set of 1000 t-tests, 50 of the tests (5 %) could be
%  false positive (<0.05). Typically, when one conducts a statistical
%  test one expects that a Type I error (one or more false positives)
%  will occur once in every 20 times the test is repeated. FWER control
%  maintains this interpretation for a set of multiple tests, such that
%  if one repeats the 1000-test experiment 20 times, at least one false
%  positive will occur in one of the experiments. FDR control permits
%  some false positives to occur at a tolerable level (e.g. 0.05). Thus,
%  each time one runs the 1000-test experiment, 5 % of the test results
%  considered significant will be declared false positive.
%
%  In the first example usage, this function applies two-tailed t-tests,
%  combined with p-value adjustments for multiplicity. The t-tests
%  compare the means of paired or unpaired measures. For unpaired
%  measures, Welch's t-test is used, which does not require the
%  condition of equal variances or equal sample sizes [1,2]. For these
%  applications, data should be given in a single vector y with groups
%  specified by a corresponding vector of group labels g (e.g. indices
%  from 1 to k). This is the general form which does not impose any
%  restriction on the number of data in each group. If y is a matrix,
%  the columns of y will be vertically concatenated and pairwise
%  comparisons will only be made between the groups within each of the
%  original y-columns but multiplicity will be adjusted across all the
%  sample groups.
%
%  This function does not support applying mixtures of paired and
%  unpaired tests. In this instance, run the data for the two groups of
%  tests seperately and create a vector of the raw p-values, as this
%  function can also be used to control the FWER or FDR for a set of
%  input p-values. In fact, because raw p-values can be used as input,
%  error rates can be controlled in almost any situation for multiple
%  testing.
%
%  The default method for controlling the FWER is the universally
%  applicable step-down Bonferroni procedure of Holm. This function
%  reports multiplicity adjusted p-values [3]. Other method options
%  available are described in detail below. They are summarised in
%  the following list, starting with the most conservative, then
%  becoming more liberal as one descends:
%
%  - Family-wise Error Rate (FWER) control:
%
%    * Holm-Bonferroni
%    * Holm-Sidak
%    * Hochberg-Bonferroni
%    * Shaffer-Bonferroni (one-way ANOVA design only)
%    * Shaffer-Sidak (one-way ANOVA design only)
%
%  - False Discovery Rate (FDR) control
%
%    * Benjamini-Hochberg
%
%
%  The input takes the following arguments
%
%   'y'       A column vector or matrix containing all the data values.
%             If fewer than two input arguments are provided, the first
%             input argument will be considered as p-values.
%
%   'g'       A column vector of the group number assignments.
%
%   'type'    The type of t-test: 'unpaired' or 'paired'. This input
%             argument is required to implement t-test comparisons.
%
%   'pairs'   A k-by-2 matrix of indices defining the comparisons,
%             where k is the number of two-sample comparisons.
%             Pairwise comparisons if set to 'all' or left empty.
%
%   'option'  The option to define the direction of the stepwise
%             Bonferroni procedure as either 'down' or 'up'. If fewer
%             than three input arguments are provided, the second input
%             argument will be considered as 'option'.
%
%             By default, the option is 'down', which implements Holm's
%             universally applicable sequentially rejective step-down
%             Bonferroni procedure [4]. In the specific circumstances
%             of complete pairwise (unpaired) comparisons across fewer
%             than 10 groups, the option is available to use Shaffer's
%             procedure, which improves on the power of Holm's method
%             by considering the logical relationships between the
%             hypotheses [5]. The algorithm uses tabulated values,
%             which are only applicable to pairwise comparisons in a
%             one-way ANOVA design [5,6]. If option is set to 'Sidak',
%             then step-down procedures are implemented, but with
%             Sidak's inequality in place of Bonferroni's. This is less
%             conservative in the typical circumstances of positive-
%             orthant dependence [6].
%
%             If the option is set to 'up', Hochberg's step-up
%             Bonferroni procedure is implemented [7]. This procedure
%             is less conservative than Holm's procedure in controlling
%             the FWER across comparisons of correlated measures.
%
%             If the option is set to 'fdr', Benjamini-Hochberg
%             procedure to control the false discovery rate is
%             implemented [8]. This procedure is less conservative
%             than the step-wise Bonferroni procedures but will not
%             maintain the FWER.
%
%
%  In the first usage example, the output is a cell array of structures
%  containing the following fields:
%
%   'test'    The tests performed.
%
%   'pair'    The pair comparison corresponding to the following
%             statistics.
%
%   'sizes'   Sample sizes.
%
%   'p'       The raw p-value of the two-sample test. These are
%             provided only for reference and should not be quoted
%             from multiple comparisons without 'padj'.
%
%   'padj'    The multiplicity adjusted p-value of the test.
%
%  In the second usage example, the output is a vector of the
%  multiplicity adjusted p-values.
%
% The algorithms used here are summarised in reference [9]
%
%
%  [1] Welch (1947) The generalization of "Student's" problem when
%       several different population variances are involved.
%       Biometrika. 34 (1???2): 28???35.
%  [2] Ruxton (2006) The unequal variance t-test is an underused
%       alternative to Student's t-test and Mann-Whitney U test.
%       Behavioral Ecology. 17(4):688.
%  [3] Wright (1992) Adjusted P-Values for Simultaneous Inference.
%       Biometrics. 48,1000-1013.
%  [4] Holm (1979) A simple sequentially rejective multiple test
%       procedure. Scandanavian Journal of Statistics, 6:65-70
%  [5] Shaffer (1986) Modified sequentially rejective multiple test
%       procedures. Journal of the American Statistical Association.
%       81: 826-831
%  [6] Holland and Copenhaver (1987) An improved sequentially
%       rejective Bonferroni test procedure. Biometrics. 43: 417-423
%  [7] Hochberg (1988) A sharper Bonferroni procedure for multiple
%       tests of significance. Biometrika. 75: 800-802
%  [8] Benjamini and Hochberg (1995) Controlling the False Discovery
%       Rate: A Practical and Powerful Approach to Multiple Testing.
%       Journal of the Royal Statistical Society. Series B. 57:289-300
%  [9] Westfall (1997) Multiple Testing of General Contrasts Using
%       Logical Constraints and Correlations. JASA 92(437):299-306
%
%  multicmp v1.0 (last updated: 21/02/2017)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%  Example data tested in Westfall (1997), JASA, 92(437): 299-306
%  p=[0.005708;0.023544;0.024193;0.044895;0.048805;0.221227;0.395867;0.693051;0.775755]

function [padj,S] = multicmp(argin1,argin2,argin3,argin4,argin5)
  if nargin<3
    % Allocate input arguments to function variables
    p = argin1;
    if size(p,2) == 1
      dim = 1;
    elseif size(p,1) == 1
      dim = 2;
    else
      error('Input p-values must be a vector');
    end
    if dim == 2
      p = p';
    end
    if nargin>1
      option = argin2;
    else
      option = 'down';
    end
    if nargout>1
      error('Invalid number of output arguments');
    end
    if size(p,1)~=numel(p)
      error('The second input argument must be a vector');
    end
    % Evaluation of test numbers for Shaffer's procedure
    k = numel(p);
    K = 3:10;
    K = K';
    n = 2+find((K.*(K-1)/2)==k);
    if strcmp(option,'down') || isempty('option')
      if ~isempty(n)
        Q = input(['Do you want to use the procedure of Shaffer for pairwise '...
                 'comparisons? (y or n) '],'s');
        if ~strcmp(Q,'y') && ~strcmp(Q,'n')
          error('The answer must be y (yes) or n (no)');
        end
      else
        Q = 'n';
      end
      if strcmp(Q,'y')
        padj = shaffer(p,'Bonferroni');
      elseif strcmp(Q,'n')
        padj = holm(p,'Bonferroni');
      end
    elseif strcmp(option,'Sidak')
      if ~isempty(n)
        Q = input(['Do you want to use the procedure of Shaffer for pairwise '...
                 'comparisons? (y or n) '],'s');
        if ~strcmp(Q,'y') && ~strcmp(Q,'n')
          error('The answer must be y (yes) or n (no)');
        end
      else
        Q = 'n';
      end
      if strcmp(Q,'y')
        padj = shaffer(p,'Sidak');
      elseif strcmp(Q,'n')
        padj = holm(p,'Sidak');
      end
    elseif strcmp(option,'up')
      padj = hochberg(p);
    elseif strcmp(option,'fdr')
      padj = fdr(p);
    end
    if dim == 2
      padj = padj';
    end
  elseif nargin > 2
    % Allocate input arguments to function variables
    y = argin1;
    g = argin2;
    type = argin3;
    if nargin>3
      pairs = argin4;
      if ~isempty(pairs)
        if strcmp(pairs,'all')
          pairs = [];
        else
          if (size(pairs,2)~=2) || ~all(all(pairs==round(pairs)))
            error('The fifth input argument must be a k-by-2 matrix of indices')
          end
        end
      end
    else
      pairs = [];
    end
    if nargin>4
      option = argin5;
    else
      option = 'down';
    end
    if nargin<3 || nargin>5
      error('Invalid number of input arguments');
    end
    if nargout>2
      error('Invalid number of output arguments');
    end
    if size(g,1)~=numel(g)
      error('The second input argument must be a column vector');
    end
    % Determine maximum number of groups defined in g-vector
    n = max(g);
    % If y is a matrix, define pairs and concatenate the y columns
    if size(y,1)~=numel(y)
      % Define pairs
      allperms = perms(1:n);
      pairs = sort(allperms(:,1:2),2);
      pairs = unique(pairs,'rows');
      pairs = pairs';
      pairs = repmat(pairs(:),1,size(y,2));
      pairs = pairs+n*repmat((0:size(y,2)-1),size(pairs,1),1);
      pairs = pairs(:);
      pairs = reshape(pairs,2,numel(pairs)/2)';
      % Redefine the vector of groups
      g = repmat(g,1,size(y,2));
      g = g+n*repmat((0:size(y,2)-1),size(y,1),1);
      % Concatenate columns of y and g
      y = y(:);
      g = g(:);
      % Redfine number of groups
      n = max(g);
    end
    % Add data to cell array
    data = cell(n,1);
    for i = 1:n
      data{i} = y(g==i);
    end
    % If not specified, determine all unique pairwise combinations
    if isempty(pairs)
      allperms = perms(1:n);
      pairs = sort(allperms(:,1:2),2);
      pairs = unique(pairs,'rows');
    end
    k = size(pairs,1);
    % Perform pairwise comparisons to obtain p-values
    p = nan(k,1);
    for i = 1:k
      if strcmp(type,'unpaired')
        p(i,1) = unpaired(data{pairs(i,1)},data{pairs(i,2)});
      elseif strcmp(type,'paired')
        p(i,1) = paired(data{pairs(i,1)},data{pairs(i,2)});
      end
    end
    % Provide strong control of family-wise error rate
    if strcmp(option,'down') || strcmp(option,'Sidak') || isempty(option)
      if (k==n*(n-1)/2) && n<=10 && n>2
        Q=input(['Do you want to use the procedure of Shaffer for pairwise '...
                 'comparisons? (y or n) '],'s');
        if ~strcmp(Q,'y') && ~strcmp(Q,'n')
          error('The answer must be y (yes) or n (no)');
        end
      else
        Q = 'n';
      end
      if strcmp(type,'unpaired') && strcmp(Q,'y')
        S = cell(k,1);
        if strcmp(option,'Sidak')
          padj = shaffer(p,'Sidak');
          for i = 1:k
            S{i}.test = 'Unpaired (Welch) t-test with Shaffer-Sidak adjustment';
          end
        else
          padj = shaffer(p,'Bonferroni');
          for i = 1:k
            S{i}.test = 'Unpaired (Welch) t-test with Shaffer-Bonferroni adjustment';
          end
        end
      else
        if strcmp(option,'Sidak')
          padj = holm(p,'Sidak');
          S = cell(k,1);
          for i = 1:k
            if strcmp(type,'unpaired')
              S{i}.test = 'Unpaired (Welch) t-test with Holm-Sidak adjustment';
            elseif strcmp(type,'paired')
              S{i}.test = 'Paired t-test with Holm-Sidak adjustment';
            end
          end
        else
          padj = holm(p,'Bonferroni');
          S = cell(k,1);
          for i=1:k
            if strcmp(type,'unpaired')
              S{i}.test = 'Unpaired (Welch) t-test with Holm-Bonferroni adjustment';
            elseif strcmp(type,'paired')
              S{i}.test = 'Paired t-test with Holm-Bonferroni adjustment';
            end
          end
        end
      end
    elseif strcmp(option,'up')
      padj = hochberg(p);
      S = cell(k,1);
      for i = 1:k
        if strcmp(type,'unpaired')
          S{i}.test = 'Unpaired (Welch) t-test with Hochberg-Bonferroni adjustment';
        elseif strcmp(type,'paired')
          S{i}.test = 'Paired t-test with Hochberg-Bonferroni adjustment';
        end
      end
    elseif strcmp(option,'fdr')
      padj = fdr(p);
      S = cell(k,1);
      for i = 1:k
        if strcmp(type,'unpaired')
          S{i}.test = 'Unpaired (Welch) t-test with Benjamini-Hochberg FDR adjustment';
        elseif strcmp(type,'paired')
          S{i}.test = 'Paired t-test with Benjamini-Hochberg FDR adjustment';
        end
      end
    end
    if nargout>1
      % Create output structure
      for i = 1:k
        S{i}.pair = pairs(i,:);
        S{i}.sizes = [numel(data{pairs(i,1)}),numel(data{pairs(i,2)})];
        S{i}.p = p(i);
        S{i}.padj = padj(i);
      end
    end
  end
end
%--------------------------------------------------------------------------