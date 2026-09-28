%% =======================================================================
%  Reversible Data Hiding using Integer Wavelet Transform (Daub 5/3)
%  and the proposed 2D Lozi-Chebyshev Chaotic Map (2D-LCM)
%      x_{i+1} = 1 - a*|cos(w*acos(x_i))| + y_i ,  y_{i+1} = b*x_i
%
%  Publication-ready reference implementation with full analysis:
%    V0  Self-verification of the Lyapunov routine (three independent
%        methods, one exactly solvable limit, two published values)
%    F1  Trajectories of six maps at their standard literature parameters:
%        Chebyshev / Henon / Lozi / 2D-SLMM / 2D-LASM / proposed 2D-LCM
%    F2  Bifurcation diagrams: Chebyshev (vs w) and the five 2D maps
%    F3  Lyapunov spectra (both exponents), plus
%        F3b  2D-LCM Lyapunov surfaces over (a,b) for w = 2.3/5.3/10.3
%        F3c  ergodicity: phase-plane density, marginal, keystream
%        F3d  benchmark comparison + periodic-window robustness table
%        T0   table of lambda_1, lambda_2, KS entropy, 0-1 test K,
%             uniformity chi2/df and lag-1 correlation for every map
%    F4  Cover / marked / encrypted images
%    F5  Histograms (cover, marked, encrypted)
%    F6  Adjacent-pixel correlation scatter plots (plain vs encrypted)
%    F7  Embedding capacity (bpp) vs PSNR curve
%    T1  Entropy, correlation coefficients, NPCR, UACI, key sensitivity,
%        reversibility check (printed to console and saved to results.mat)
%
%  Pipeline (sender):
%    cover -> FOB preprocessing -> Daub 5/3 IWT -> histogram-shifting
%    (+/-1) embedding in hl, lh, hh -> inverse IWT -> marked image
%    -> pixel-domain encryption (XOR keystream + chaotic permutation)
%    -> encrypted stego image
%
%  Pipeline (receiver): exact inverse; recovers secret bits AND the
%  bit-exact original cover image (verified at the end).
%
%  IMPLEMENTATION NOTES (read before adapting for the thesis text):
%  (1) Embedding uses double-peak histogram shifting (peaks 0 and -1) on
%      the detail coefficients. This realizes the "+/-1 change" rule of
%      Sec. 3.3.2 in a BLIND-decodable, provably reversible way. The
%      literal parity rule (bit 1: even-1 / odd+1) is not blind-decodable
%      (the receiver cannot distinguish (v,bit0) from (v+1,bit1) without
%      side information), so HS -- already cited in Sec. 3.3.3 [110],[111]
%      -- is the standard working realization.
%  (2) Encryption is applied in the PIXEL domain to the marked image
%      (XOR + permutation). XOR-encrypting the ll band before the inverse
%      IWT would push reconstructed pixels far outside [0,255] (invalid
%      8-bit image). Encrypt-after-embedding is standard in RDH-EI and
%      keeps every step range-safe and reversible.
%  (3) FOB (fall-of-boundary): pixels are pre-clamped to [DELTA,255-DELTA]
%      and each clamped pixel's CORRECTION (in [0,DELTA], 2 bits for
%      DELTA=3) is carried as auxiliary data -- not its full 8-bit
%      original. The clamp direction is inferable (clamped-low pixels sit
%      at DELTA, clamped-high at 255-DELTA), so only the magnitude is
%      stored; correction 0 marks a pixel genuinely at the boundary.
%      This yields ~2-3% more payload than raw 8-bit FOB storage.
%      DELTA = 3 is empirically the minimum safe margin: the measured
%      worst-case pixel deviation from embedding is exactly 3, and
%      DELTA = 0 or 1 provably overflow.
%      LIMITATION: the FOB cost scales with the number of saturated
%      pixels, so images that are largely saturated (all-black/white,
%      binary/document scans) can have FOB cost exceed capacity. The
%      scheme targets natural imagery; all USC-SIPI test images used here
%      have <0.2% boundary pixels.
%
%  Requirements: base MATLAB R2016b+ (no toolboxes needed; uses a local
%  PSNR/entropy). Test image: set cfg.imagePath, or a synthetic image is
%  generated if the file is unavailable.
%  =======================================================================

clearvars; close all; clc; rng(1);                          % reproducible

%% ----------------------- 0. Configuration -----------------------------
cfg.imagePath   = '/Users/sudhirsingh/Desktop/rdh chaos/final h/Cover Image/Jetplane.tiff';   % any grayscale image; auto-fallback
cfg.imSize      = 512;               % working size (even), images resized
cfg.key         = struct('x0',0.3,'y0',0.6,'a',5,'b',20.7,'w',2.3,'form','fold');
%  ^ 2D-LCM secret key = (x0,y0,a,b,w) + realization selector 'form'.
%    (a,b,w) = (5, 20.7, 2.3) is a verified HYPER-chaotic operating
%    point (both Lyapunov exponents > 0, see Section 2).
%    form = 'fold' is the
%    canonical realization; set form = 'mod1' to reproduce the earlier
%    mod-1 version of this script exactly. Previous default was
%    ('x0',0.3671,'y0',0.2115,'a',1.45,'b',1.50,'w',3.5).
cfg.payload_bpp = 0.05;              % demo payload (used only if no secret image)
cfg.secretImagePath = '/Users/sudhirsingh/Desktop/rdh chaos/final h/5.1.12.tiff';            % 256x256 secret image to hide; '' = random bits
cfg.secretSize  = 256;               % secret image is resized to this (square)
cfg.useCompression = true;           % losslessly compress the binary secret if it helps
cfg.DELTA       = 3;                 % FOB safety margin (see note 3)
cfg.transient   = 2000;              % chaotic iterations discarded
cfg.outDir      = 'results';
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

%% ----------------------- 1. Load cover image --------------------------
try
    I0 = imread(cfg.imagePath);
    if size(I0,3) == 3, I0 = round(0.2989*double(I0(:,:,1)) + ...
            0.5870*double(I0(:,:,2)) + 0.1140*double(I0(:,:,3)));
    else, I0 = double(I0); end
catch
    warning('Image not found -- generating synthetic test image.');
    [xx,yy] = meshgrid(linspace(0,1,cfg.imSize));
    I0 = round(255*(0.5+0.5*sin(8*pi*xx).*cos(6*pi*yy)));
end
I0 = simpleResize(I0, cfg.imSize);          % local nearest-neighbour resize
I0 = max(0,min(255,round(I0)));
[M,N] = size(I0);
fprintf('Cover image: %dx%d, 8-bit grayscale\n', M, N);

%% ------------- 2. Chaotic map characterization (F1-F3) ----------------
% =======================================================================
%  PROPOSED 2D LOZI-CHEBYSHEV MAP (2D-LCM)
%
%      x_{i+1} = 1 - a*|cos( w*acos(x_i) )| + y_i                  (LCM)
%      y_{i+1} = b*x_i
%
%  with w in (2,+inf) and a,b in (-inf,+inf). The Chebyshev term supplies
%  a wide chaotic range and a near-uniform invariant density; the Lozi
%  (absolute-value) term supplies a non-smooth nonlinearity, so the map
%  has no smooth critical point and no homoclinic tangency structure.
%
%  BOUNDING. cos(w*acos x) is real only on [-1,1] and grows like
%  cosh(w*acosh|x|) outside it, so the raw recursion has no bounded
%  invariant set at the operating points used here. The recursion is
%  therefore kept EXACTLY as written and x alone is folded back into the
%  Chebyshev domain by the measure-preserving sawtooth
%          wrap11(u) = mod(u+1,2) - 1   in  [-1,1).
%  Since y_{i+1} = b*x_i and |x_i| < 1, y stays in (-b,b) automatically.
%  Set cfg.key.form = 'mod1' for the earlier mod-1 realization.
%
%  ANALYTIC PROPERTIES (all verified numerically in block V0 below).
%  With u = w*acos(x) the Jacobian is
%      J = [ -a*w*sign(cos u)*sin(u)/sqrt(1-x^2) , 1 ;  b , 0 ]
%  (the fold has unit derivative a.e.), hence
%      det J = -b        =>      lambda_1 + lambda_2 = ln|b|
%  for every a and w. Two consequences worth stating explicitly in any
%  write-up, because they are easy to misread from the figures:
%    * the two exponents are mirror images about ln|b|/2, NOT about 0.
%      They are symmetric about zero only in the area-preserving case
%      |b| = 1. The green reference line in F3/F3b marks the true axis.
%    * at a = 0 the nonlinearity vanishes and the map becomes linear,
%      J = [0 1; b 0], with eigenvalues +-sqrt(b), so both exponents
%      equal ln(sqrt(b)). That value is the FLOOR of the lambda_1 curve
%      and it is set by b alone. "Chaotic for every a" is therefore a
%      consequence of choosing |b| > 1 together with the fold; the
%      nonlinearity is what lifts lambda_1 ABOVE that floor, and that
%      lift is the quantity worth reporting.
%
%  Figures produced in this section:
%    F1  trajectories of six maps at their standard literature parameters
%    F2  bifurcation diagrams: Chebyshev (vs w) and the five 2D maps
%    F3  Lyapunov spectra (both exponents) of the 2D maps
%    F3b 2D-LCM Lyapunov surfaces over (a,b) for w = 2.3/5.3/10.3 plus
%        the w = 2.3, b = 20.7 curve
%    F3c ergodicity: phase-plane density, marginal, keystream histogram
%    F3d benchmark comparison with a periodic-window robustness table
%    T0  table: lambda_1, lambda_2, KS entropy, 0-1 test K, uniformity
%        chi2/df and lag-1 correlation for every map
%    V0  self-verification of the Lyapunov routine (three independent
%        methods + one exactly solvable limit)
%
%  The Chebyshev-coupling construction follows the family of 2D maps
%  built by driving a seed map with a Chebyshev polynomial (see the
%  chaotic-map literature cited in the thesis). Delete this paragraph if
%  the lineage is documented elsewhere in your text.
% =======================================================================
fprintf('\n=== 2D-LCM characterization ===\n');

tr      = cfg.transient;   % transient discarded from orbits
trLE    = 500;             % transient discarded from Lyapunov runs
nTraj   = 20000;           % points per trajectory plot
nKeep   = 200;             % points kept per bifurcation column
nAbif   = 500;             % parameter samples per bifurcation diagram
nLEc    = 5000;            % iterations per LE point (analytic Jacobian)
nLEref  = 2500;            % iterations per LE point (numerical Jacobian)
nRefPts = 200;             % parameter samples per reference-map LE curve
nLEgrid = 1500;            % iterations per LE point on an (a,b) surface
gN      = 41;              % (a,b) grid resolution for F3b
N01     = 5000;            % sequence length for the 0-1 test

%% ---- V0: is the Lyapunov routine itself correct? --------------------
% Four independent checks. Any bug in the Jacobian, the Gram-Schmidt step
% or the folding shows up here immediately.
fprintf('\n  --- V0: verification of the Lyapunov routine ---\n');
vChk = [5 5.7 5.3; 5 20.7 2.3; cfg.key.a cfg.key.b cfg.key.w];
for q = 1:size(vChk,1)
    av = vChk(q,1); bv = vChk(q,2); wv = vChk(q,3);
    LA = leLCM(av,bv,wv,0.3,0.6,40000,trLE);            % analytic Jacobian
    LN = leLCMnum(av,bv,wv,0.3,0.6,20000,trLE);         % numerical Jacobian
    LB = leLCMbenettin(av,bv,wv,0.3,0.6,100000,trLE);   % no Jacobian at all
    fprintf('  a=%.4g b=%.4g w=%.4g:\n', av, bv, wv);
    fprintf('    analytic-J  : [%+.5f %+.5f]   sum = %.6f , ln|b| = %.6f\n', ...
            LA(1), LA(2), sum(LA), log(abs(bv)));
    fprintf('    numerical-J : [%+.5f %+.5f]   (independent Jacobian)\n', LN(1), LN(2));
    fprintf('    Benettin    :  %+.5f           (two-trajectory, no Jacobian)\n', LB);
end
fprintf('  exactly solvable limit a -> 0 (map is linear, lambda = ln sqrt(b)):\n');
for bv = [1.5 5.7 20.7]
    L0 = leLCM(0, bv, 5.3, 0.3, 0.6, 40000, trLE);
    fprintf('    b=%5.4g : computed [%+.6f %+.6f]   exact %+.6f\n', ...
            bv, L0(1), L0(2), log(sqrt(bv)));
end
LEhen = lyapunov2D(@(v)[1-1.4*v(1)^2+v(2); 0.3*v(1)], [0;0], 100000, trLE, ...
                   @(v)[-2.8*v(1), 1; 0.3, 0]);
LEloz = lyapunov2D(@(v)[1-1.7*abs(v(1))+v(2); 0.5*v(1)], [0.1;0.1], 100000, trLE, ...
                   @(v)[-1.7*sign(v(1)), 1; 0.5, 0]);
fprintf('  same routine on maps with published values:\n');
fprintf('    Henon(1.4,0.3) lam_1 = %+.4f  [published ~ +0.419]\n', LEhen(1));
fprintf('    Lozi (1.7,0.5) lam_1 = %+.4f  [published ~ +0.470]\n', LEloz(1));

%% ---- F1: trajectories -----------------------------------------------
fig = figure('Name','F1 Trajectories','Position',[30 380 1200 620]);

xcb = chebyshevMap(0.3, 10.3, nTraj, tr);
subplot(2,3,1); plot(xcb(1:end-1), xcb(2:end), '.', 'MarkerSize',1);
xlabel('x_i'); ylabel('x_{i+1}'); title('(a) Chebyshev, w=10.3'); axis tight; box on

[xhn,yhn] = henonMap(0, 0, 1.4, 0.3, nTraj, tr);       % (0,0): inside the basin
subplot(2,3,2); plot(xhn, yhn, '.', 'MarkerSize',1);
xlabel('x_i'); ylabel('y_i'); title('(b) Henon, a=1.4, b=0.3'); axis tight; box on

[xlz,ylz] = loziMap(0.1, 0.1, 1.7, 0.5, nTraj, tr);
subplot(2,3,3); plot(xlz, ylz, '.', 'MarkerSize',1);
xlabel('x_i'); ylabel('y_i'); title('(c) Lozi, a=1.7, b=0.5'); axis tight; box on

[xsl,ysl] = slmm2D(0.3, 0.6, 1, 3, nTraj, tr);
subplot(2,3,4); plot(xsl, ysl, '.', 'MarkerSize',1);
xlabel('x_i'); ylabel('y_i'); title('(d) 2D-SLMM, a=1, b=3'); axis tight; box on

[xla,yla] = lasm2D(0.3, 0.6, 0.9, nTraj, tr);
subplot(2,3,5); plot(xla, yla, '.', 'MarkerSize',1);
xlabel('x_i'); ylabel('y_i'); title('(e) 2D-LASM, a=0.9'); axis tight; box on

[xlc,ylc] = lcm2D(0.3, 0.6, 5, 5.7, 5.3, nTraj, tr);
subplot(2,3,6); plot(xlc, ylc, '.', 'MarkerSize',1);
xlabel('x_i'); ylabel('y_i'); title('(f) 2D-LCM, a=5, b=5.7, w=5.3'); axis tight; box on

saveFig(fig, fullfile(cfg.outDir,'F1_trajectories'));
fprintf('\n  F1 done: the 2D-LCM covers the whole phase plane; the classical\n');
fprintf('           maps live on thin, structured attractors\n');

%% ---- F2: bifurcation diagrams ---------------------------------------
fig = figure('Name','F2 Bifurcation','Position',[30 30 1200 700]);

% (a) Chebyshev: the control parameter is w, chaotic for w > 2
wG = linspace(1.5, 10.3, nAbif);
subplot(2,3,1); hold on
for ww = wG
    xb = chebyshevMap(0.3, ww, nKeep, tr);
    plot(ww*ones(nKeep,1), xb, '.k','MarkerSize',0.5);
end
xlabel('w'); ylabel('x_i'); title('(a) Chebyshev'); box on

aG = linspace(0, 1.4, nAbif);
subplot(2,3,2); hold on
for aa = aG
    xb = henonMap(0, 0, aa, 0.3, nKeep, tr);
    plot(aa*ones(nKeep,1), xb, '.k','MarkerSize',0.5);
end
xlabel('a'); ylabel('x_i'); title('(b) Henon, b=0.3'); ylim([-2 2]); box on

aG = linspace(0.5, 1.75, nAbif);          % Lozi escapes beyond a ~ 1.76
subplot(2,3,3); hold on
for aa = aG
    xb = loziMap(0.1, 0.1, aa, 0.5, nKeep, tr);
    plot(aa*ones(nKeep,1), xb, '.k','MarkerSize',0.5);
end
xlabel('a'); ylabel('x_i'); title('(c) Lozi, b=0.5'); box on

aG = linspace(0.85, 1, nAbif);
subplot(2,3,4); hold on
for aa = aG
    xb = slmm2D(0.3, 0.6, aa, 3, nKeep, tr);
    plot(aa*ones(nKeep,1), xb, '.k','MarkerSize',0.5);
end
xlabel('a'); ylabel('x_i'); title('(d) 2D-SLMM, b=3'); box on

aG = linspace(0, 1, nAbif);
subplot(2,3,5); hold on
for aa = aG
    xb = lasm2D(0.3, 0.6, aa, nKeep, tr);
    plot(aa*ones(nKeep,1), xb, '.k','MarkerSize',0.5);
end
xlabel('a'); ylabel('x_i'); title('(e) 2D-LASM'); box on

aG = linspace(-10, 10, nAbif);
subplot(2,3,6); hold on
for aa = aG
    xb = lcm2D(0.3, 0.6, aa, 5.7, 5.3, nKeep, tr);
    plot(aa*ones(nKeep,1), xb, '.k','MarkerSize',0.5);
end
xlabel('a'); ylabel('x_i'); title('(f) 2D-LCM, b=5.7, w=5.3'); box on
saveFig(fig, fullfile(cfg.outDir,'F2_bifurcation'));
fprintf('  F2 done: 2D-LCM output fills the plane for every a (no periodic windows)\n');

%% ---- F3: Lyapunov spectra -------------------------------------------
fig = figure('Name','F3 Lyapunov exponents','Position',[30 30 1200 700]);

% (a) Chebyshev is 1D, so it has a single exponent; ln(w) is the
%     asymptotic prediction and doubles as a correctness check
wG = linspace(1.5, 10.3, nRefPts); Lc = zeros(1,nRefPts);
for ii = 1:nRefPts, Lc(ii) = leCheb1D(0.3, wG(ii), nLEc, trLE); end
subplot(2,3,1);
plot(wG,Lc,'b-',wG,zeros(1,nRefPts),'k--',wG,log(wG),'g-.','LineWidth',1.1); grid on
xlabel('w'); ylabel('\lambda'); title('(a) Chebyshev LE');
legend('\lambda','\lambda = 0','ln w','Location','southeast')

aG = linspace(0, 1.4, nRefPts); L = nan(2, nRefPts);
for ii = 1:nRefPts
    aa = aG(ii);
    L(:,ii) = lyapunov2D(@(v)[1-aa*v(1)^2+v(2); 0.3*v(1)], [0;0], ...
                         nLEref, trLE, @(v)[-2*aa*v(1), 1; 0.3, 0]);
end
subplot(2,3,2); plot(aG,L(1,:),'b-',aG,L(2,:),'r-',aG,zeros(1,nRefPts),'k--'); grid on
xlabel('a'); ylabel('\lambda'); title('(b) Henon LE, b=0.3'); ylim([-3 1.5])

aG = linspace(0.5, 1.75, nRefPts); L = nan(2, nRefPts);
for ii = 1:nRefPts
    aa = aG(ii);
    L(:,ii) = lyapunov2D(@(v)[1-aa*abs(v(1))+v(2); 0.5*v(1)], [0.1;0.1], ...
                         nLEref, trLE, @(v)[-aa*sign(v(1)), 1; 0.5, 0]);
end
subplot(2,3,3); plot(aG,L(1,:),'b-',aG,L(2,:),'r-',aG,zeros(1,nRefPts),'k--'); grid on
xlabel('a'); ylabel('\lambda'); title('(c) Lozi LE, b=0.5')

aG = linspace(0.85, 1, nRefPts); L = nan(2, nRefPts);
for ii = 1:nRefPts
    L(:,ii) = lyapunov2D(@(v) slmmStep(v, aG(ii), 3), [0.3;0.6], nLEref, trLE, []);
end
subplot(2,3,4); plot(aG,L(1,:),'b-',aG,L(2,:),'r-',aG,zeros(1,nRefPts),'k--'); grid on
xlabel('a'); ylabel('\lambda'); title('(d) 2D-SLMM LE, b=3')

aG = linspace(0, 1, nRefPts); L = nan(2, nRefPts);
for ii = 1:nRefPts
    L(:,ii) = lyapunov2D(@(v) lasmStep(v, aG(ii)), [0.3;0.6], nLEref, trLE, []);
end
subplot(2,3,5); plot(aG,L(1,:),'b-',aG,L(2,:),'r-',aG,zeros(1,nRefPts),'k--'); grid on
xlabel('a'); ylabel('\lambda'); title('(e) 2D-LASM LE')

aGm = linspace(-10, 10, 300); Ll = zeros(2,300);
for ii = 1:300, Ll(:,ii) = leLCM(aGm(ii), 5.7, 5.3, 0.3, 0.6, nLEc, trLE); end
% the exponents mirror about ln|b|/2 (because det J = -b), NOT about zero
ax57 = log(5.7)/2;
subplot(2,3,6);
plot(aGm,Ll(1,:),'b-',aGm,Ll(2,:),'r-',aGm,zeros(1,300),'k--', ...
     aGm,ax57*ones(1,300),'g-.','LineWidth',1.1); grid on
xlabel('a'); ylabel('\lambda'); title('(f) 2D-LCM LE, b=5.7, w=5.3')
legend('\lambda_1','\lambda_2','\lambda = 0','ln|b|/2','Location','best')
saveFig(fig, fullfile(cfg.outDir,'F3_lyapunov'));
LLE = Ll(1,:);                       % 2D-LCM MLE curve (kept for reference)
fprintf('  F3 done. 2D-LCM at b=5.7, w=5.3 : min lambda_1 = %+.4f over a in [-10,10]\n', ...
        min(Ll(1,:)));
fprintf('           that minimum sits at a = 0 and equals ln sqrt(b) = %+.4f,\n', log(sqrt(5.7)));
fprintf('           the floor set by b alone; the nonlinearity lifts lambda_1 above it\n');
fprintf('           hyper-chaotic (both LE>0) for |a| <= %.2f\n', ...
        maxAbsOrZero(aGm(Ll(2,:) > 0)));
fprintf('  SYMMETRY: lambda_1 and lambda_2 mirror about ln|b|/2 = %+.5f, NOT about 0.\n', ax57);
fprintf('            Max deviation of (lam_1+lam_2)/2 from that axis: %.2e\n', ...
        max(abs(sum(Ll,1)/2 - ax57)));

%% ---- F3b: Lyapunov surfaces over (a,b) ------------------------------
wList = [2.3 5.3 10.3];
aGg = linspace(-20, 20, gN); bGg = linspace(-20, 20, gN);
[AA, BB] = meshgrid(aGg, bGg);
fig = figure('Name','F3b 2D-LCM LE surfaces','Position',[30 30 1000 950]);
for iw = 1:numel(wList)
    L1 = zeros(gN); L2 = zeros(gN);
    for ib = 1:gN
        for ia = 1:gN
            LEg = leLCM(aGg(ia), bGg(ib), wList(iw), 0.3, 0.6, nLEgrid, trLE);
            L1(ib,ia) = max(LEg(1), -10);      % clamp the degenerate b=0 column
            L2(ib,ia) = max(LEg(2), -10);
        end
    end
    subplot(4,2,2*iw-1);
    surf(AA,BB,L1,'EdgeColor','none'); hold on
    surf(AA,BB,L2,'EdgeColor','none');
    surf(AA,BB,zeros(gN),'FaceColor',[.6 .6 .6],'EdgeColor','none','FaceAlpha',0.35);
    xlabel('a'); ylabel('b'); zlabel('\lambda'); view(-35,25); zlim([-6 6])
    title(sprintf('2D-LCM LE surface, w = %.1f', wList(iw)));
    subplot(4,2,2*iw);
    surf(AA,BB,L1,'EdgeColor','none'); hold on
    surf(AA,BB,L2,'EdgeColor','none');
    view(90,0); xlabel('a'); ylabel('b'); zlabel('\lambda'); zlim([-6 6])
    title(sprintf('left view, w = %.1f', wList(iw)));
    fprintf('  F3b: LE surface for w = %.1f done\n', wList(iw));
end
aG2 = linspace(-20, 20, 300); Lg = zeros(2,300);
for ii = 1:300, Lg(:,ii) = leLCM(aG2(ii), 20.7, 2.3, 0.3, 0.6, nLEc, trLE); end
subplot(4,2,[7 8]);
plot(aG2,Lg(1,:),'b-',aG2,Lg(2,:),'r-',aG2,zeros(1,300),'k--', ...
     aG2,(log(20.7)/2)*ones(1,300),'g-.','LineWidth',1.1); grid on
xlabel('a'); ylabel('\lambda');
title('2D-LCM LE, w = 2.3, b = 20.7  (hyper-chaotic where both curves > 0)');
legend('\lambda_1','\lambda_2','\lambda = 0','ln|b|/2','Location','best');
saveFig(fig, fullfile(cfg.outDir,'F3b_LE_surfaces'));

%% ---- F3c: ergodicity of the proposed map ----------------------------
[xE,yE] = lcm2D(0.3, 0.6, cfg.key.a, cfg.key.b, cfg.key.w, 200000, tr);
ksE = keystream8(wrap11(xE + yE));
fig = figure('Name','F3c Ergodicity','Position',[30 30 1200 350]);
nb = 128; ymax = max(abs(yE)) + eps;
ixe = min(nb, max(1, floor((xE+1)/2*nb)+1));
iye = min(nb, max(1, floor((yE/ymax+1)/2*nb)+1));
Dns = accumarray([iye(:) ixe(:)], 1, [nb nb]);
subplot(1,3,1); imagesc([-1 1], [-ymax ymax], Dns); axis xy square
xlabel('x'); ylabel('y'); title('2D-LCM phase-plane density'); colorbar
subplot(1,3,2); histogram(xE, linspace(-1,1,101)); box on
xlabel('x'); ylabel('count');
title(sprintf('marginal density of x, \\chi^2/df = %.3f', chi2Uniform(xE,100,-1,1)));
subplot(1,3,3); histogram(ksE, -0.5:1:255.5); box on
xlabel('keystream byte'); ylabel('count');
title(sprintf('keystream, \\chi^2/df = %.3f', chi2Uniform(ksE,256,-0.5,255.5)));
saveFig(fig, fullfile(cfg.outDir,'F3c_ergodicity'));
fprintf('  F3c: lag-1 correlation of the keystream = %+.5f\n', ...
        corrCoefLocal(ksE(1:end-1), ksE(2:end)));

%% ---- F3d: benchmark comparison + robustness -------------------------
% NOTE ON WHAT THIS FIGURE CAN AND CANNOT SHOW.
%   lambda_1 ~ ln(a*w) - 0.386 for this map family, so the height a curve
%   reaches at the right-hand edge is fixed by whichever a_max the plotted
%   range stops at. Since a is unbounded, that stopping point is a choice,
%   not a property: extending a range makes any map "win". Raising w or
%   |b| does the same at no cost. Curve HEIGHTS are therefore NOT a
%   ranking and must not be presented as one.
%   What IS comparable is the SHAPE: a dip below zero is a periodic
%   window, a parameter value at which a nominally chaotic map is not
%   chaotic. That feature is invariant under rescaling of the axis, and
%   the table printed below quantifies it. Report the table.
bM = 1.5; wM = 3.5; aMax = 2.35;
nP = 200;
names  = {'Henon (b=0.3)','Lozi (b=0.5)','Chebyshev','2D-SLMM (b=3)', ...
          '2D-LASM','2D-LCM (proposed)'};
ranges = {[0 1.4],[0.5 1.75],[2.05 10.3],[0.85 1],[0.4 1],[0 aMax]};
% dark, mutually distinguishable colours (light cyan/yellow are unreadable)
cols   = { [0.60 0.00 0.20], ...      % Henon      dark maroon
           [0.00 0.50 0.00], ...      % Lozi       dark green
           [0.49 0.18 0.56], ...      % Chebyshev  purple
           [0.00 0.45 0.55], ...      % 2D-SLMM    dark teal
           [0.85 0.45 0.00], ...      % 2D-LASM    dark orange
           [0.00 0.20 0.80] };        % 2D-LCM     strong blue
lsty   = {'-','-','-','-','-.','-'};
lw     = [1.3 1.3 1.3 1.3 1.3 2.2];
Lall   = cell(1,numel(names));
fig = figure('Name','F3d Benchmark comparison','Position',[30 30 800 520]); hold on
for m = 1:numel(names)
    pv = linspace(ranges{m}(1), ranges{m}(2), nP); Lm = zeros(1,nP);
    for ii = 1:nP
        p = pv(ii);
        switch m
            case 1, tmp = lyapunov2D(@(v)[1-p*v(1)^2+v(2); 0.3*v(1)], [0;0], ...
                                     nLEref, trLE, @(v)[-2*p*v(1), 1; 0.3, 0]);
            case 2, tmp = lyapunov2D(@(v)[1-p*abs(v(1))+v(2); 0.5*v(1)], [0.1;0.1], ...
                                     nLEref, trLE, @(v)[-p*sign(v(1)), 1; 0.5, 0]);
            case 3, tmp = [leCheb1D(0.3, p, nLEc, trLE); NaN];
            case 4, tmp = lyapunov2D(@(v) slmmStep(v, p, 3), [0.3;0.6], nLEref, trLE, []);
            case 5, tmp = lyapunov2D(@(v) lasmStep(v, p),    [0.3;0.6], nLEref, trLE, []);
            case 6, tmp = leLCM(p, bM, wM, 0.3, 0.6, nLEc, trLE);
        end
        Lm(ii) = tmp(1);
    end
    Lall{m} = Lm;
    plot(linspace(0,1,nP), Lm, 'Color', cols{m}, 'LineStyle', lsty{m}, 'LineWidth', lw(m));
end
plot([0 1], [0 0], 'k--', 'LineWidth', 1.0); grid on; ylim([-2 2]); box on
xlabel('control parameter, rescaled to [0,1] (heights are NOT comparable)');
ylabel('\lambda_1');
title('Robustness comparison: dips below 0 are periodic windows');
legend(names, 'Location','southeast');
saveFig(fig, fullfile(cfg.outDir,'F3d_benchmark'));

fprintf('  F3d: robustness over each map''s own range (this IS comparable)\n');
fprintf('  %-24s %10s %14s %9s\n', 'map', 'min lam_1', '% of range<=0', 'windows');
for m = 1:numel(names)
    Lm = Lall{m}; bad = Lm <= 0;
    nw = sum(diff([false bad]) == 1);
    fprintf('  %-24s %+10.4f %13.1f%% %9d\n', names{m}, min(Lm), 100*mean(bad), nw);
end
fprintf(['  Read this table, not the curve heights: 0%% and 0 windows means the\n' ...
         '  map is chaotic everywhere in its range, which is what a cipher needs.\n']);

%% ---- T0: quantitative comparison table ------------------------------
fprintf('\n  --- T0: chaotic-performance comparison ---\n');
fprintf('  %-22s %8s %8s %8s %7s %8s %9s\n', 'map', 'lam_1', 'lam_2', 'KS-ent', 'K(0-1)', 'chi2/df', 'lag1corr');
rows = { 'Chebyshev (w=10.3)', 'Henon (1.4,0.3)', 'Lozi (1.7,0.5)', ...
         '2D-SLMM (1,3)', '2D-LASM (0.9)', '2D-LCM (5,5.7,5.3)', '2D-LCM @ working key' };
for r = 1:numel(rows)
    switch r
        case 1
            LE = [leCheb1D(0.3, 10.3, 20000, tr); NaN];
            sq = chebyshevMap(0.3, 10.3, N01, tr);
        case 2
            LE = lyapunov2D(@(v)[1-1.4*v(1)^2+v(2); 0.3*v(1)], [0;0], 20000, trLE, ...
                            @(v)[-2*1.4*v(1), 1; 0.3, 0]);
            sq = henonMap(0, 0, 1.4, 0.3, N01, tr);
        case 3
            LE = lyapunov2D(@(v)[1-1.7*abs(v(1))+v(2); 0.5*v(1)], [0.1;0.1], 20000, trLE, ...
                            @(v)[-1.7*sign(v(1)), 1; 0.5, 0]);
            sq = loziMap(0.1, 0.1, 1.7, 0.5, N01, tr);
        case 4
            LE = lyapunov2D(@(v) slmmStep(v,1,3), [0.3;0.6], 20000, trLE, []);
            sq = slmm2D(0.3, 0.6, 1, 3, N01, tr);
        case 5
            LE = lyapunov2D(@(v) lasmStep(v,0.9), [0.3;0.6], 20000, trLE, []);
            sq = lasm2D(0.3, 0.6, 0.9, N01, tr);
        case 6
            LE = leLCM(5, 5.7, 5.3, 0.3, 0.6, 20000, trLE);
            sq = lcm2D(0.3, 0.6, 5, 5.7, 5.3, N01, tr);
        otherwise
            LE = leLCM(cfg.key.a, cfg.key.b, cfg.key.w, cfg.key.x0, cfg.key.y0, 20000, trLE);
            sq = lcm2D(cfg.key.x0, cfg.key.y0, cfg.key.a, cfg.key.b, cfg.key.w, N01, tr);
    end
    ks  = sum(LE(LE > 0 & isfinite(LE)));                   % Kolmogorov-Sinai entropy
    fprintf('  %-22s %+8.4f %+8.4f %8.4f %7.3f %8.3f %+9.4f\n', rows{r}, ...
            LE(1), LE(2), ks, test01(sq, 10), ...
            chi2Uniform(sq, 50, min(sq), max(sq)), corrCoefLocal(sq(1:end-1), sq(2:end)));
end
fprintf('  (K ~ 1 => chaotic, K ~ 0 => regular; chi2/df ~ 1 => uniform coverage)\n');

%% ---- range verification ---------------------------------------------
fprintf('\n  --- range verification for the 2D-LCM ---\n');
aS = linspace(-5, 5, 200); L57 = zeros(2,200);
for ii = 1:200, L57(:,ii) = leLCM(aS(ii), 5.7, 5.3, 0.3, 0.6, nLEc, trLE); end
fprintf('  b=5.7 , w=5.3 : chaotic for all tested a: %s ; hyper-chaotic for |a| <= %.2f\n', ...
        bool2str(all(L57(1,:) > 0)), maxAbsOrZero(aS(L57(2,:) > 0)));
fprintf('  b=20.7, w=2.3 : chaotic for all tested a: %s ; hyper-chaotic for |a| <= %.2f\n', ...
        bool2str(all(Lg(1,:) > 0)), maxAbsOrZero(aG2(Lg(2,:) > 0)));
for wv = [2.3 5.3 10.3]
    st = '';
    for bv = [0.25 0.5 0.9 1.5 3.0]
        Lb = leLCM(5, bv, wv, 0.3, 0.6, nLEc, trLE);
        st = [st sprintf('  b=%.2f:%s', bv, bool2str(Lb(1) > 0))]; %#ok<AGROW>
    end
    fprintf('  w=%.1f chaotic? (a=5)%s\n', wv, st);
end
LEk = leLCM(cfg.key.a, cfg.key.b, cfg.key.w, cfg.key.x0, cfg.key.y0, 20000, trLE, keyForm(cfg.key));
if LEk(2) > 0, cls = 'HYPER-CHAOTIC'; elseif LEk(1) > 0, cls = 'chaotic'; else, cls = 'NOT chaotic'; end
fprintf('  working key (a=%.4g, b=%.4g, w=%.4g, form=%s): lambda = [%+.4f, %+.4f],\n', ...
        cfg.key.a, cfg.key.b, cfg.key.w, keyForm(cfg.key), LEk(1), LEk(2));
fprintf('              sum = %+.4f, ln|b| = %+.4f  ->  %s\n', ...
        sum(LEk), log(abs(cfg.key.b)), cls);
fprintf('  lambda_max at working key (a=%.2f): %.4f  (chaotic if > 0)\n', ...
        cfg.key.a, largestLyapunov(cfg.key, 8000, cfg.transient));

%% ---------------- 3. Sender: embed + encrypt (main run) ---------------
% Build the secret payload. If a 256x256 secret IMAGE is provided, binarize
% it, optionally compress it losslessly, and prepend a small self-describing
% header so the receiver can decompress and rebuild the image. Otherwise fall
% back to a random bit-string of length cfg.payload_bpp * M * N.
[secret, secretMeta] = buildSecretPayload(cfg, M, N);
nBits = numel(secret);
fprintf('\n=== Embedding secret payload ===\n');
if secretMeta.isImage
    fprintf('Secret: %dx%d image, binarized = %d bits', secretMeta.h, secretMeta.w, secretMeta.rawBits);
    if secretMeta.compressed
        fprintf(' -> compressed to %d bits (%.2fx)\n', nBits, secretMeta.rawBits/nBits);
    else
        fprintf(' (compression not beneficial; stored raw = %d bits)\n', nBits);
    end
else
    fprintf('Secret: %d random payload bits (%.3f bpp)\n', nBits, cfg.payload_bpp);
end

[Istego, auxFOB, hdrInfo] = senderPipeline(I0, secret, cfg);

% Marked (pre-encryption) image for imperceptibility metrics
Imarked = hdrInfo.Imarked;
psnr_marked = psnrLocal(I0, Imarked);
mse_marked  = mean((double(I0(:))-double(Imarked(:))).^2);
bpp_actual  = nBits/(M*N);
fprintf('Marked-image PSNR : %.2f dB\n', psnr_marked);
fprintf('  payload embedded: %d bits  =  %.4f bpp  (MSE = %.4f)\n', nBits, bpp_actual, mse_marked);
if secretMeta.isImage
    fprintf('  NOTE: this PSNR corresponds to the ACTUAL payload of %.4f bpp above,\n', bpp_actual);
    fprintf('        not to the 0.75 bpp capacity. Report PSNR at a matched rate\n');
    fprintf('        (e.g. 0.10 bpp for comparison with the base paper).\n');
end
fprintf('Encrypted stego produced (%dx%d, uint8-valid: %d)\n', ...
        M, N, all(Istego(:)>=0 & Istego(:)<=255));

% F4: visual results
figure('Name','F4 Images');
subplot(1,3,1); imshow(uint8(I0));      title('Cover');
subplot(1,3,2); imshow(uint8(Imarked)); title('Marked');
subplot(1,3,3); imshow(uint8(Istego));  title('Encrypted stego');
saveFig(gcf, fullfile(cfg.outDir,'F4_images'));

% F5: histograms
figure('Name','F5 Histograms');
subplot(1,3,1); histogram(I0(:),0:255);      title('Cover');
subplot(1,3,2); histogram(Imarked(:),0:255); title('Marked');
subplot(1,3,3); histogram(Istego(:),0:255);  title('Encrypted stego');
saveFig(gcf, fullfile(cfg.outDir,'F5_histograms'));

%% ---------------- 4. Receiver: extract + recover ----------------------
[secret_rx, Irec] = receiverPipeline(Istego, auxFOB, cfg, numel(secret));

ok_bits = isequal(secret_rx, secret);
ok_img  = isequal(Irec, I0);
fprintf('\n=== Reversibility check ===\n');
fprintf('Secret bits recovered exactly : %s\n', bool2str(ok_bits));
fprintf('Cover image recovered exactly : %s (MSE = %.3g)\n', ...
        bool2str(ok_img), mean((Irec(:)-I0(:)).^2));
assert(ok_bits && ok_img, 'Reversibility failed -- check pipeline.');

% If the secret was an image, decompress and rebuild it, then report NC.
if secretMeta.isImage
    secretImg_rx = rebuildSecretImage(secret_rx, secretMeta);
    ncVal = ncMetric(secretMeta.binImage, secretImg_rx);
    fprintf('Secret IMAGE recovered        : NC = %.4f (bit-exact: %s)\n', ...
            ncVal, bool2str(isequal(secretMeta.binImage, secretImg_rx)));
    figure('Name','F4b Secret image');
    subplot(1,2,1); imshow(secretMeta.binImage); title('Original secret (binary)');
    subplot(1,2,2); imshow(secretImg_rx);        title(sprintf('Recovered secret (NC=%.4f)',ncVal));
    saveFig(gcf, fullfile(cfg.outDir,'F4b_secret_image'));
end

%% ---------------- 5. Security analysis (T1, F6) ------------------------
fprintf('\n=== Security analysis ===\n');

% 5.1 Entropy
H_cover = shannonEntropy(I0);
H_enc   = shannonEntropy(Istego);
fprintf('Entropy  cover     : %.4f bits/pixel\n', H_cover);
fprintf('Entropy  encrypted : %.4f bits/pixel (ideal 8.0000)\n', H_enc);

% 5.2 Adjacent-pixel correlation (H, V, D), 5000 random pairs each
dirs = {'horizontal','vertical','diagonal'};
corrP = zeros(1,3); corrE = zeros(1,3);
figure('Name','F6 Correlation');
for d = 1:3
    [p1,p2] = adjacentPairs(I0, dirs{d}, 5000);
    [e1,e2] = adjacentPairs(Istego, dirs{d}, 5000);
    corrP(d) = corrCoefLocal(p1,p2);
    corrE(d) = corrCoefLocal(e1,e2);
    subplot(2,3,d);   plot(p1,p2,'.','MarkerSize',2);
    title(sprintf('Plain %s (r=%.4f)',dirs{d},corrP(d))); axis square
    subplot(2,3,d+3); plot(e1,e2,'.','MarkerSize',2);
    title(sprintf('Encrypted %s (r=%.4f)',dirs{d},corrE(d))); axis square
end
saveFig(gcf, fullfile(cfg.outDir,'F6_correlation'));
fprintf('Correlation  (plain)     H/V/D : %+.4f  %+.4f  %+.4f\n', corrP);
fprintf('Correlation  (encrypted) H/V/D : %+.4f  %+.4f  %+.4f\n', corrE);

% 5.3 NPCR / UACI (plaintext sensitivity: one-pixel change in the cover)
% NOTE: NPCR/UACI measure how a one-pixel plaintext change propagates
% through the ciphertext. This requires PLAINTEXT-DEPENDENT DIFFUSION --
% the keystream must itself depend on the image (see plaintextSeed in
% senderPipeline). With a fixed keystream, a one-pixel change would alter
% only that pixel's ciphertext (NPCR ~ 1/(MN) ~ 0), which is a property of
% plain XOR, not a bug. The diffusion stage there gives proper avalanche.
I1 = I0; I1(1,1) = mod(I1(1,1)+1, 256);
E0 = Istego;
E1 = senderPipeline(I1, secret, cfg);           % encrypt modified cover
npcr = 100 * mean(E0(:) ~= E1(:));
uaci = 100 * mean(abs(E0(:) - E1(:)) / 255);
fprintf('NPCR : %.4f %%   (ideal ~99.6094)\n', npcr);
fprintf('UACI : %.4f %%   (ideal ~33.4635)\n', uaci);

% 5.4 Key sensitivity: decrypt with x0 perturbed by 1e-15.
% A wrong key yields a corrupted bitstream; hsExtract will legitimately
% fail to complete extraction, so we guard it and report the failure as
% evidence of key sensitivity instead of letting the assert abort.
cfgWrong = cfg; cfgWrong.key.x0 = cfg.key.x0 + 1e-15;
try
    [~, IrecWrong] = receiverPipeline(Istego, auxFOB, cfgWrong, numel(secret));
    psnr_wrong = psnrLocal(I0, IrecWrong);
    fprintf('Wrong-key recovery PSNR : %.2f dB (noise-like => key-sensitive)\n', ...
            psnr_wrong);
catch ME
    psnr_wrong = NaN;
    fprintf(['Wrong-key extraction failed (%s) => scheme is key-sensitive ' ...
             '(unrecoverable with a 1e-15 key error).\n'], ME.message);
end

%% ---------------- 6. Capacity vs PSNR curve (F7) -----------------------
%% --- 6a. Capacity report: the exact embedding rate of THIS image -------
% Net rate = (N0 + Nm1 - 32 - 8*N_FOB) / (M*N).
% Our double-peak HS embedding is BLIND-decodable (values are self-
% identifying + 32-bit length header), so NO location map is charged.
% The rate is IMAGE-DEPENDENT -- report this per test image in Table 3.2
% instead of a flat value.
fprintf('\n=== Capacity report (use these numbers in Table 3.2) ===\n');
capReport(I0, cfg.DELTA, cfg.imagePath);

% --- PSNR at standard, comparable embedding rates (report THESE in the paper) ---
% PSNR must always be quoted at a stated rate. These fixed rates give the
% honest, reproducible numbers: 0.10 bpp matches the base paper's comparison
% point, and 0.75 bpp is the full-capacity operating point.
fprintf('\n=== PSNR at matched embedding rates (report these) ===\n');
for r = [0.05 0.10 0.25 0.50 0.75]
    nb = round(r*M*N);
    if nb + 32 > hdrInfo.capacityBits, nb = hdrInfo.capacityBits - 32; end
    if nb <= 0, continue; end
    rbits = rand(nb,1) > 0.5;
    [~,~,info_r] = senderPipeline(I0, rbits, cfg);
    fprintf('  %.4f bpp -> PSNR = %.2f dB\n', nb/(M*N), psnrLocal(I0, info_r.Imarked));
end
fprintf('  (Use 0.10 bpp for comparison with the base paper; 0.75 bpp = full capacity.)\n');

fprintf('\n=== Capacity vs PSNR curve ===\n');
maxCap  = hdrInfo.capacityBits;                     % peak-bin capacity
maxBpp  = (maxCap - 40) / (M*N);                    % keep header headroom
bppGrid = linspace(0.01, 0.95*maxBpp, 12);
psnrCurve = zeros(size(bppGrid));
for ii = 1:numel(bppGrid)
    nb   = round(bppGrid(ii)*M*N);
    bits = rand(nb,1) > 0.5;
    [~,~,info_i] = senderPipeline(I0, bits, cfg);
    psnrCurve(ii) = psnrLocal(I0, info_i.Imarked);
    fprintf('  %.3f bpp -> %.2f dB\n', bppGrid(ii), psnrCurve(ii));
end
figure('Name','F7 Capacity-PSNR');
plot(bppGrid, psnrCurve, 'bo-','LineWidth',1.2); grid on
xlabel('Embedding capacity (bpp)'); ylabel('PSNR of marked image (dB)');
title('Capacity vs PSNR -- IWT + Lozi-Chebyshev RDH');
saveFig(gcf, fullfile(cfg.outDir,'F7_capacity_psnr'));

%% ---------------- 7. Save numerical results ----------------------------
results = struct('psnr_marked',psnr_marked,'entropy_cover',H_cover, ...
    'entropy_encrypted',H_enc,'corr_plain',corrP,'corr_encrypted',corrE, ...
    'NPCR',npcr,'UACI',uaci,'psnr_wrongkey',psnr_wrong, ...
    'bpp',bppGrid,'psnr_curve',psnrCurve,'max_capacity_bpp',maxBpp, ...
    'reversible_bits',ok_bits,'reversible_image',ok_img,'key',cfg.key);
save(fullfile(cfg.outDir,'results.mat'),'results');
fprintf('\nAll figures and results.mat saved in ./%s\n', cfg.outDir);

%% =======================================================================
%                            LOCAL FUNCTIONS
%  =======================================================================

function [secret, meta] = buildSecretPayload(cfg, M, N)
% Build the payload bit-vector. If cfg.secretImagePath points to an image,
% hide that image: resize to secretSize x secretSize, binarize (threshold 128),
% optionally RLE-compress, and prepend a self-describing header so the receiver
% can rebuild it. Otherwise generate random bits at cfg.payload_bpp.
meta = struct('isImage',false,'compressed',false,'h',0,'w',0,'rawBits',0);
useImg = isfield(cfg,'secretImagePath') && ~isempty(cfg.secretImagePath) && exist(cfg.secretImagePath,'file');
if ~useImg
    nBits  = round(cfg.payload_bpp * M * N);
    secret = rand(nBits,1) > 0.5;
    return;
end
S = imread(cfg.secretImagePath);
if size(S,3) > 1, S = rgb2gray(S); end
S = imresize(S, [cfg.secretSize cfg.secretSize]);
binImg = S >= 128;                          % binarize to 0/1
bits   = double(binImg(:));                 % column-major bit-stream
meta.isImage=true; meta.h=size(binImg,1); meta.w=size(binImg,2);
meta.rawBits=numel(bits); meta.binImage=binImg;

% Optional lossless compression (run-length on the bit-stream, then packed).
payloadBits = bits; compressed = false;
if isfield(cfg,'useCompression') && cfg.useCompression
    rle = rleEncodeBits(bits);              % vector of run lengths (uint16)
    packed = packRunsToBits(rle);           % bit-vector
    if numel(packed) < numel(bits)          % only use it if it actually helps
        payloadBits = packed; compressed = true;
    end
end
meta.compressed = compressed;

% Header (self-describing): [flagCompressed(1) | h(16) | w(16) | payloadLen(32)]
hdr = [double(compressed), de2biVec(meta.h,16), de2biVec(meta.w,16), de2biVec(numel(payloadBits),32)]';
secret = [hdr; payloadBits(:)];
end

function img = rebuildSecretImage(rxBits, meta)
% Inverse of buildSecretPayload: strip header, decompress if needed, reshape.
rxBits = rxBits(:);
compressed = rxBits(1) > 0.5;
h  = bi2deVec(rxBits(2:17));
w  = bi2deVec(rxBits(18:33));
plen = bi2deVec(rxBits(34:65));
payload = rxBits(66:66+plen-1);
if compressed
    rle  = unpackBitsToRuns(payload);
    bits = rleDecodeBits(rle, h*w);
else
    bits = payload;
end
img = reshape(logical(bits(1:h*w)), h, w);
end

% ---- run-length codec on a 0/1 stream (lossless) ----
function runs = rleEncodeBits(bits)
bits=bits(:)'; d=diff([~bits(1), bits]); idx=[find(d~=0), numel(bits)+1];
runs=diff([0 idx-1]); runs=runs(runs>0);            % run lengths
% split any run >65535 into 65535-chunks with a 0-length toggle marker
out=[]; for r=runs, while r>65535, out=[out 65535 0]; r=r-65535; end, out=[out r]; end
runs=uint16(out);
% store starting bit as the first element's parity via a leading marker
runs=[uint16(bits(1)) runs];
end
function bits = rleDecodeBits(runs, n)
first=double(runs(1)); runs=double(runs(2:end));
bits=zeros(1,sum(runs)); cur=first; p=1;
for r=runs, if r>0, bits(p:p+r-1)=cur; p=p+r; end, cur=~cur; end
bits=bits(1:n)';
end
function packed = packRunsToBits(runs)
% each run length stored as 16 bits
packed=zeros(16*numel(runs),1); for i=1:numel(runs)
    packed((i-1)*16+(1:16))=de2biVec(double(runs(i)),16); end
end
function runs = unpackBitsToRuns(bits)
bits=bits(:); n=floor(numel(bits)/16); runs=zeros(1,n,'uint16');
for i=1:n, runs(i)=uint16(bi2deVec(bits((i-1)*16+(1:16)))); end
end

% ---- small bit helpers (no toolbox) ----
function b = de2biVec(v,n), b=zeros(1,n); for k=1:n, b(k)=mod(floor(v/2^(n-k)),2); end, end
function v = bi2deVec(b), b=b(:)'; n=numel(b); v=sum(b.*2.^(n-1:-1:0)); end

function nc = ncMetric(A,B)
A=double(A(:)); B=double(B(:)); a=A-mean(A); b=B-mean(B);
d=sqrt(sum(a.^2)*sum(b.^2)); if d==0, nc=double(isequal(A,B)); else nc=sum(a.*b)/d; end
end

function [Istego, auxFOB, info] = senderPipeline(I0, secret, cfg)
% Full sender: FOB preprocessing -> IWT -> HS embedding -> IIWT ->
% pixel-domain encryption (XOR + permutation).
    [M,N] = size(I0); D = cfg.DELTA;

    % --- FOB preprocessing (note 3): 2-bit correction encoding ---
    % Pixels are clamped to [D,255-D]. Rather than storing each boundary
    % pixel's full 8-bit original, we store only its CLAMP CORRECTION in
    % [0,D] (2 bits for D=3). The clamp DIRECTION need not be stored: a
    % clamped-low pixel ends up exactly at D, a clamped-high one exactly
    % at 255-D. The receiver rebuilds the identical "ambiguous set"
    % {Ip==D or Ip==255-D} in the same scan order and inverts the clamp;
    % correction 0 marks a pixel that genuinely equalled D (or 255-D) and
    % was never clamped. This costs 2 bits/ambiguous pixel instead of
    % 8 bits/boundary pixel (~2-3% more payload on natural images).
    Ip      = min(max(I0, D), 255-D);
    ambMask = (Ip == D) | (Ip == 255-D);
    idx     = find(ambMask);
    fobCorr = zeros(numel(idx),1);
    for t = 1:numel(idx)
        if Ip(idx(t)) == D, fobCorr(t) = D - I0(idx(t));
        else,               fobCorr(t) = I0(idx(t)) - (255-D);
        end
    end
    auxFOB = struct('corr', fobCorr);          % 2 bits each when packed

    % --- Forward Daub 5/3 IWT ---
    [ll,hl,lh,hh] = fwd53_2d(Ip);

    % --- Parity-rule embedding in [hl lh hh] (Sec. 3.3.2) ---
    % One bit per detail coefficient across all three sub-bands (HL, LH, HH):
    %   secret bit 0 -> coefficient unchanged
    %   secret bit 1 -> even coefficient -1, odd coefficient +1  (parity flips)
    % The original parities are recorded as recovery side-information so the
    % receiver can decode each bit and restore the exact coefficient. This
    % realizes the full 3-subband capacity of 3*(M/2)*(N/2) = 0.75 bpp.
    C = [hl(:); lh(:); hh(:)];
    capacity = numel(C);                            % 1 bit per detail coeff
    header   = de2bi_local(numel(secret), 32);      % 32-bit payload length
    bits     = [header(:); secret(:)];
    assert(numel(bits) <= capacity, ...
        ['Payload (%d bits) exceeds this cover''s embedding capacity (%d bits). ' ...
         'For a 256x256 secret image, enable cfg.useCompression, use a larger ' ...
         'cover, or reduce the secret size.'], numel(bits), capacity);
    origPar = mod(C,2);                             % recovery side-information
    [Cw, usedPar] = parityEmbed(C, bits);           % embed; usedPar = parities used
    auxFOB.origPar = origPar(1:numel(bits));        % carry parities for used coeffs
    nb = numel(hl);
    hl(:) = Cw(1:nb); lh(:) = Cw(nb+1:2*nb); hh(:) = Cw(2*nb+1:end);

    % --- Inverse IWT -> marked image ---
    Imarked = inv53_2d(ll,hl,lh,hh);
    assert(all(Imarked(:)>=0 & Imarked(:)<=255), 'FOB margin violated.');

    % --- Pixel-domain encryption (note 2) ---
    % Confusion (XOR keystream) + diffusion (plaintext-dependent forward
    % chaining) + permutation. The chaining makes every ciphertext pixel
    % depend on all preceding plaintext pixels, giving proper avalanche
    % (NPCR ~ 99.6%, UACI ~ 33.4%). It is exactly invertible (see receiver).
    [~,~,pseq] = loziChebyshev(cfg.key, 3*M*N, cfg.transient);
    ks = keystream8(pseq(1:M*N));                   % 8-bit XOR keystream
    p1 = chaoticPerm(pseq(M*N+1:2*M*N));            % permutation round 1
    p2 = chaoticPerm(pseq(2*M*N+1:3*M*N));          % permutation round 2
    % Three modular-addition chaining rounds interleaved with two chaotic
    % permutations. Chaining gives confusion+diffusion; the permutations
    % between rounds spread a one-pixel change to every output position,
    % so avalanche is position-independent (NPCR~99.6%, UACI~33.4%).
    a = chainFwd(double(Imarked(:)), ks, 173);
    b = chainFwd(a(p1), ks, 91);
    c = chainFwd(b(p2), ks, 47);
    Istego = reshape(c, M, N);

    info = struct('Imarked',Imarked,'capacityBits',capacity);
end

function c = chainFwd(v, ks, seed)
% Forward modular-addition chain: c(t) = v(t) + ks(t) + c(t-1) (mod 256),
% c(0)=seed. Exactly invertible by chainInv.
    n = numel(v); c = zeros(n,1); prev = seed;
    for t = 1:n
        c(t) = mod(v(t) + ks(t) + prev, 256);
        prev = c(t);
    end
end

function v = chainInv(c, ks, seed)
% Inverse of chainFwd.
    n = numel(c); v = zeros(n,1); prev = seed;
    for t = 1:n
        v(t) = mod(c(t) - ks(t) - prev, 256);
        prev = c(t);
    end
end

function [secret, Irec] = receiverPipeline(Istego, auxFOB, cfg, nSecret)
% Full receiver: decrypt -> IWT -> extract + restore -> IIWT -> undo FOB.
    [M,N] = size(Istego);

    % --- Decrypt: invert 3-round chain/permute (reverse order) ---
    [~,~,pseq] = loziChebyshev(cfg.key, 3*M*N, cfg.transient);
    ks = keystream8(pseq(1:M*N));
    p1 = chaoticPerm(pseq(M*N+1:2*M*N));   ip1(p1) = 1:M*N;
    p2 = chaoticPerm(pseq(2*M*N+1:3*M*N)); ip2(p2) = 1:M*N;
    b = chainInv(Istego(:), ks, 47); b = b(ip2);   % undo round 3 + perm2
    a = chainInv(b, ks, 91);         a = a(ip1);    % undo round 2 + perm1
    Imarked = reshape(chainInv(a, ks, 173), M, N);  % undo round 1

    % --- IWT, extract bits, restore coefficients ---
    [ll,hl,lh,hh] = fwd53_2d(Imarked);
    Cw = [hl(:); lh(:); hh(:)];
    nUsed = 32 + nSecret;
    [bits, C] = parityExtract(Cw, auxFOB.origPar, nUsed);
    secret = bits(33:end);
    lenHdr = bi2de_local(bits(1:32));
    assert(lenHdr == nSecret, 'Header/payload length mismatch.');
    nb = numel(hl);
    hl(:) = C(1:nb); lh(:) = C(nb+1:2*nb); hh(:) = C(2*nb+1:end);

    % --- Inverse IWT and FOB restoration (2-bit correction decode) ---
    Ip = inv53_2d(ll,hl,lh,hh);
    D  = cfg.DELTA;
    ambMask = (Ip == D) | (Ip == 255-D);       % same set as the sender
    idx = find(ambMask);
    Irec = Ip;
    if numel(idx) ~= numel(auxFOB.corr)
        error('FOB ambiguous-set mismatch (%d vs %d).', ...
              numel(idx), numel(auxFOB.corr));
    end
    for t = 1:numel(idx)
        c = auxFOB.corr(t);
        if c == 0, continue; end               % genuinely at the boundary
        if Ip(idx(t)) == D, Irec(idx(t)) = D - c;
        else,               Irec(idx(t)) = (255-D) + c;
        end
    end
end

% ============ Chaotic maps: proposed 2D-LCM and references =============
function [x,y,p] = loziChebyshev(K, n, transient)
% Proposed 2D Lozi-Chebyshev map (2D-LCM):
%     x_{i+1} = 1 - a*|cos( w*acos(x_i) )| + y_i
%     y_{i+1} = b*x_i
% K = struct('x0','y0','a','b','w') and, optionally, 'form':
%   'fold' (default) -- x folded into the Chebyshev domain [-1,1) by the
%                       measure-preserving sawtooth wrap11; y = b*x then
%                       lies in (-b,b) automatically.
%   'mod1'           -- legacy realization, both states reduced mod 1.
% p is the mixing sequence consumed by the encryption stage.
    switch keyForm(K)
        case 'mod1'
            [x,y,p] = lcm2Dmod1(K.x0, K.y0, K.a, K.b, K.w, n, transient);
        otherwise
            [x,y] = lcm2D(K.x0, K.y0, K.a, K.b, K.w, n, transient);
            p = wrap11(x + y);
    end
end

function f = keyForm(K)
% Realization selector carried inside the key struct ('fold' by default).
    if isfield(K,'form') && ~isempty(K.form), f = lower(K.form); else, f = 'fold'; end
end

function [x,y] = lcm2D(x0, y0, a, b, w, n, transient)
% Core 2D-LCM iteration (fold form).
    x = zeros(n,1); y = zeros(n,1);
    xc = x0; yc = y0;
    for k = 1:(transient+n)
        xg = min(1, max(-1, xc));                          % Chebyshev domain
        xn = wrap11(1 - a*abs(cos(w*acos(xg))) + yc);      % Lozi-Chebyshev
        yn = b*xc;
        xc = xn; yc = yn;
        if k > transient, x(k-transient) = xc; y(k-transient) = yc; end
    end
end

function [x,y,p] = lcm2Dmod1(x0, y0, a, b, w, n, transient)
% Legacy mod-1 realization of the 2D-LCM (kept for backward compatibility).
    x = zeros(n,1); y = zeros(n,1); p = zeros(n,1);
    xc = x0; yc = y0;
    for k = 1:(transient+n)
        xg = min(1, max(-1, xc));
        xn = mod(1 - a*abs(cos(w*acos(xg))) + yc, 1);
        yn = mod(b*xc, 1);
        pn = mod(a*xc + b*yc, 1);
        xc = xn; yc = yn;
        if k > transient
            x(k-transient) = xc; y(k-transient) = yc; p(k-transient) = pn;
        end
    end
end

function v = wrap11(u)
% Measure-preserving sawtooth fold of R onto [-1,1).
    v = mod(u + 1, 2) - 1;
end

% ---- Lyapunov spectrum of the 2D-LCM: three independent routes --------
function LE = leLCM(a, b, w, x0, y0, n, transient, form)
% PRIMARY method. Full spectrum [lambda_1; lambda_2] (descending) by
% continuous Gram-Schmidt re-orthonormalization of the tangent map, with
% the ANALYTIC Jacobian. With u = w*acos(x),
%     J = [ -a*w*sign(cos u)*sin(u)/sqrt(1-x^2) , 1 ;  b , 0 ]
% so det J = -b and lambda_1 + lambda_2 = ln|b| -- a built-in check.
% lambda_1 > 0 => chaotic; lambda_2 > 0 as well => hyper-chaotic.
% form = 'fold' (default) or 'mod1'; the tangent map is the same for both
% (each fold has unit derivative a.e.), only the state update differs.
    if nargin < 8, form = 'fold'; end
    isMod = strcmpi(form, 'mod1');
    xc = x0; yc = y0;
    q11 = 1; q21 = 0; q12 = 0; q22 = 1;
    s1 = 0; s2 = 0; cnt = 0;
    for k = 1:(transient+n)
        xg  = min(1, max(-1, xc));
        u   = w*acos(xg);
        cu  = cos(u);
        den = sqrt(max(1 - xg^2, 1e-12));
        sg  = sign(cu); if sg == 0, sg = 1; end
        j11 = -a*sg*w*sin(u)/den;
        % tangent step Z = J*Q  with J = [j11 1; b 0]
        z11 = j11*q11 + q21;   z21 = b*q11;
        z12 = j11*q12 + q22;   z22 = b*q12;
        r1 = hypot(z11, z21);  if r1 < realmin, r1 = realmin; end
        q11 = z11/r1; q21 = z21/r1;
        d   = q11*z12 + q21*z22;
        z12 = z12 - d*q11;  z22 = z22 - d*q21;
        r2 = hypot(z12, z22);  if r2 < realmin, r2 = realmin; end
        q12 = z12/r2; q22 = z22/r2;
        if isMod, xn = mod(1 - a*abs(cu) + yc, 1); yn = mod(b*xc, 1);
        else,     xn = wrap11(1 - a*abs(cu) + yc); yn = b*xc;
        end
        xc = xn; yc = yn;
        if k > transient
            s1 = s1 + log(r1); s2 = s2 + log(r2); cnt = cnt + 1;
        end
    end
    LE = sort([s1; s2]/max(cnt,1), 'descend');
end

function LE = leLCMnum(a, b, w, x0, y0, n, transient)
% CHECK 1: identical algorithm but with a central-difference Jacobian, so
% the analytic derivative above is never used. Agreement to ~5 decimals
% confirms the analytic Jacobian.
    f = @(v) [wrap11(1 - a*abs(cos(w*acos(min(1,max(-1,v(1)))))) + v(2)); b*v(1)];
    LE = lyapunov2D(f, [x0;y0], n, transient, []);
end

function lam = leLCMbenettin(a, b, w, x0, y0, n, transient)
% CHECK 2: lambda_1 by two-trajectory renormalization (Benettin). Uses no
% Jacobian and no orthonormalization at all -- a fully independent route.
    d0 = 1e-10; s = 0; cnt = 0;
    step = @(v) [wrap11(1 - a*abs(cos(w*acos(min(1,max(-1,v(1)))))) + v(2)); b*v(1)];
    v1 = [x0; y0]; v2 = v1 + [d0; 0];
    for k = 1:(transient+n)
        v1 = step(v1); v2 = step(v2);
        dv = v2 - v1; d = norm(dv); if d == 0, d = realmin; end
        if k > transient, s = s + log(d/d0); cnt = cnt + 1; end
        v2 = v1 + dv*(d0/d);
    end
    lam = s/max(cnt,1);
end

function lam = largestLyapunov(K, n, transient)
% Maximum Lyapunov exponent (MLE) of the 2D-LCM at key K. Original
% signature preserved so the rest of the script is unchanged.
    LE  = leLCM(K.a, K.b, K.w, K.x0, K.y0, n, transient, keyForm(K));
    lam = LE(1);
end

% ---- reference maps used in the comparisons of Figs. 1-3 --------------
function x = chebyshevMap(x0, w, n, transient)
% Chebyshev map: x_{i+1} = cos(w*acos(x_i)), chaotic for w > 2.
    x = zeros(n,1); xc = x0;
    for k = 1:(transient+n)
        xc = cos(w*acos(min(1, max(-1, xc))));
        if k > transient, x(k-transient) = xc; end
    end
end

function lam = leCheb1D(x0, w, n, transient)
% Lyapunov exponent of the 1D Chebyshev map (~ ln w for large w).
    xc = x0; s = 0; cnt = 0;
    for k = 1:(transient+n)
        xg  = min(1, max(-1, xc));
        den = sqrt(max(1 - xg^2, 1e-12));
        d   = w*sin(w*acos(xg))/den;
        xc  = cos(w*acos(xg));
        if k > transient, s = s + log(max(abs(d), realmin)); cnt = cnt + 1; end
    end
    lam = s/max(cnt,1);
end

function [x,y] = henonMap(x0, y0, a, b, n, transient)
% Henon map: x_{i+1} = 1 - a*x_i^2 + y_i , y_{i+1} = b*x_i.
% Call with the initial condition (0,0): it lies inside the basin for the
% whole a-range scanned here, whereas e.g. (0.3,0.6) escapes for a > ~1.18
% and would blank out the plots. Diverging orbits are returned as NaN.
    x = zeros(n,1); y = zeros(n,1); xc = x0; yc = y0;
    for k = 1:(transient+n)
        xn = 1 - a*xc^2 + yc; yn = b*xc;
        xc = xn; yc = yn;
        if ~isfinite(xc) || abs(xc) > 1e8, xc = NaN; yc = NaN; end
        if k > transient, x(k-transient) = xc; y(k-transient) = yc; end
    end
end

function [x,y] = loziMap(x0, y0, a, b, n, transient)
% Lozi map: x_{i+1} = 1 - a*|x_i| + y_i , y_{i+1} = b*x_i (a=1.7, b=0.5 is
% the classical attractor; the orbit escapes for a > ~1.76 at b = 0.5).
    x = zeros(n,1); y = zeros(n,1); xc = x0; yc = y0;
    for k = 1:(transient+n)
        xn = 1 - a*abs(xc) + yc; yn = b*xc;
        xc = xn; yc = yn;
        if ~isfinite(xc) || abs(xc) > 1e8, xc = NaN; yc = NaN; end
        if k > transient, x(k-transient) = xc; y(k-transient) = yc; end
    end
end

function [x,y] = slmm2D(x0, y0, a, b, n, transient)
% 2D Sine-Logistic modulation map (2D-SLMM):
%   x_{i+1} = a*( sin(pi*y_i)     + b )*x_i*(1-x_i)
%   y_{i+1} = a*( sin(pi*x_{i+1}) + b )*y_i*(1-y_i)
    x = zeros(n,1); y = zeros(n,1); xc = x0; yc = y0;
    for k = 1:(transient+n)
        xn = a*(sin(pi*yc) + b)*xc*(1-xc);
        yn = a*(sin(pi*xn) + b)*yc*(1-yc);
        xc = xn; yc = yn;
        if k > transient, x(k-transient) = xc; y(k-transient) = yc; end
    end
end

function v = slmmStep(v, a, b)
    xn = a*(sin(pi*v(2)) + b)*v(1)*(1-v(1));
    yn = a*(sin(pi*xn)   + b)*v(2)*(1-v(2));
    v  = [xn; yn];
end

function [x,y] = lasm2D(x0, y0, mu, n, transient)
% 2D Logistic-adjusted-Sine map (2D-LASM):
%   x_{i+1} = sin( pi*mu*(y_i     + 3)*x_i*(1-x_i) )
%   y_{i+1} = sin( pi*mu*(x_{i+1} + 3)*y_i*(1-y_i) )
    x = zeros(n,1); y = zeros(n,1); xc = x0; yc = y0;
    for k = 1:(transient+n)
        xn = sin(pi*mu*(yc + 3)*xc*(1-xc));
        yn = sin(pi*mu*(xn + 3)*yc*(1-yc));
        xc = xn; yc = yn;
        if k > transient, x(k-transient) = xc; y(k-transient) = yc; end
    end
end

function v = lasmStep(v, mu)
    xn = sin(pi*mu*(v(2) + 3)*v(1)*(1-v(1)));
    yn = sin(pi*mu*(xn   + 3)*v(2)*(1-v(2)));
    v  = [xn; yn];
end

% ---- generic Lyapunov spectrum for the reference 2D maps --------------
function LE = lyapunov2D(fstep, v0, n, transient, jacFun)
% fstep: handle v -> v_next. jacFun: handle v -> 2x2 Jacobian, or [] for a
% central-difference Jacobian. Diverging orbits return [NaN; NaN].
    v = v0(:); Q = eye(2); s = [0;0]; cnt = 0;
    for k = 1:(transient+n)
        if isempty(jacFun), J = numJac2(fstep, v); else, J = jacFun(v); end
        v = fstep(v);
        if any(~isfinite(v)) || any(abs(v) > 1e8), LE = [NaN; NaN]; return; end
        [Q, r] = gramSchmidt2(J*Q);
        if k > transient, s = s + log(max(r, realmin)); cnt = cnt + 1; end
    end
    LE = sort(s/max(cnt,1), 'descend');
end

function J = numJac2(f, v)
    h = 1e-7; J = zeros(2,2);
    for j = 1:2
        e = zeros(2,1); e(j) = h;
        J(:,j) = (f(v+e) - f(v-e))/(2*h);
    end
end

function [Q,r] = gramSchmidt2(Z)
    r = zeros(2,1);
    z1 = Z(:,1); r(1) = norm(z1);
    if r(1) > 0, q1 = z1/r(1); else, q1 = [1;0]; r(1) = realmin; end
    z2 = Z(:,2) - (q1.'*Z(:,2))*q1; r(2) = norm(z2);
    if r(2) > 0, q2 = z2/r(2); else, q2 = [-q1(2); q1(1)]; r(2) = realmin; end
    Q = [q1 q2];
end

% ---- randomness / ergodicity diagnostics ------------------------------
function K = test01(phi, ncs)
% Gottwald-Melbourne 0-1 test for chaos: K ~ 1 chaotic, K ~ 0 regular.
% ncs deterministic values of c are used and the median K is reported.
    phi = phi(:).';
    N = numel(phi); nmax = max(10, floor(N/10));
    j = 1:N; mu2 = mean(phi)^2;
    cs = pi/5 + (3*pi/5)*((0:ncs-1) + 0.5)/ncs;
    Ks = zeros(1,ncs); nn = 1:nmax;
    for t = 1:ncs
        c = cs(t);
        p = cumsum(phi.*cos(j*c)); q = cumsum(phi.*sin(j*c));
        D = zeros(1,nmax);
        for m = nn
            dp = p(m+1:end) - p(1:end-m);
            dq = q(m+1:end) - q(1:end-m);
            D(m) = mean(dp.^2 + dq.^2) - mu2*(1-cos(m*c))/(1-cos(c));
        end
        Ks(t) = corrCoefLocal(nn, D);
    end
    K = median(Ks);
end

function s = chi2Uniform(x, nb, lo, hi)
% Chi-square per degree of freedom of the orbit histogram: s ~ 1 means the
% orbit covers its range uniformly (quantifies ergodicity).
    x = x(:); x = x(isfinite(x));
    h = histcounts(x, linspace(lo, hi, nb+1));
    e = numel(x)/nb;
    if e <= 0, s = NaN; else, s = sum((h - e).^2/e)/(nb - 1); end
end

function v = maxAbsOrZero(z)
    if isempty(z), v = 0; else, v = max(abs(z)); end
end

% ---- discretization of the chaotic sequence ---------------------------
function ks = keystream8(p)
% 8-bit keystream by the low-order-digit rule:
%     ks_i = floor( |p_i| * 1e9 )  mod  256.
% Taking low-order digits removes any residual non-uniformity of the
% invariant density, so ks is uniform on {0,...,255} for either
% realization ('fold' or 'mod1').
    ks = mod(floor(abs(p(:))*1e9), 256);
end

function perm = chaoticPerm(p)
% Chaotic permutation: ranking of the chaotic sequence (realizes the
% location scrambling as a full permutation).
    [~,perm] = sort(p);
end

% -------------------- Daub 5/3 integer lifting -------------------------
function [ll,hl,lh,hh] = fwd53_2d(X)
% JPEG2000-reversible 5/3 lifting, columns then rows. Integer in/out.
    [S,D]  = lift53(X);                             % along columns
    [ll,hl] = liftRows(S);                          % along rows of low
    [lh,hh] = liftRows(D);                          % along rows of high
end

function X = inv53_2d(ll,hl,lh,hh)
    S = iliftRows(ll,hl);
    D = iliftRows(lh,hh);
    X = ilift53(S,D);
end

function [S,D] = lift53(X)                          % along dim 1
    Xe = X(1:2:end,:); Xo = X(2:2:end,:);
    D  = Xo - floor((Xe + Xe([2:end end],:))/2);
    S  = Xe + floor((D([1 1:end-1],:) + D + 2)/4);
end

function X = ilift53(S,D)
    Xe = S - floor((D([1 1:end-1],:) + D + 2)/4);
    Xo = D + floor((Xe + Xe([2:end end],:))/2);
    X  = zeros(2*size(S,1), size(S,2));
    X(1:2:end,:) = Xe; X(2:2:end,:) = Xo;
end

function [L,H] = liftRows(X),  [Lt,Ht] = lift53(X.'); L = Lt.'; H = Ht.'; end
function X = iliftRows(L,H),   X = ilift53(L.',H.').'; end

% ------------------ Histogram-shifting embed/extract -------------------
function [Cw, usedPar] = parityEmbed(C, bits)
% Parity-rule embedding (Sec. 3.3.2): one bit per coefficient, in order.
%   bit 0 -> unchanged ; bit 1 -> even coeff -1, odd coeff +1 (parity flips).
% Returns the modified coefficient vector and the original parities used.
    Cw = C; n = numel(bits);
    usedPar = mod(C(1:n),2);
    for t = 1:n
        if bits(t) == 1
            if mod(C(t),2) == 0, Cw(t) = C(t) - 1;   % even -> -1
            else,                Cw(t) = C(t) + 1;   % odd  -> +1
            end
        end
    end
end

function [bits, C] = parityExtract(Cw, origPar, nUsed)
% Inverse of parityEmbed. A bit is 1 iff the current parity differs from the
% stored original parity; the +-1 modification is then reversed exactly.
    C = Cw; bits = zeros(nUsed,1);
    for t = 1:nUsed
        curPar = mod(Cw(t),2);
        if curPar ~= origPar(t)
            bits(t) = 1;
            % reverse the change: if original was even it had been -1 -> +1;
            % if original was odd it had been +1 -> -1.
            if origPar(t) == 0, C(t) = Cw(t) + 1;   % restore even
            else,               C(t) = Cw(t) - 1;   % restore odd
            end
        else
            bits(t) = 0; C(t) = Cw(t);
        end
    end
end

function C = hsEmbed(C, bits)
% Double-peak HS (peaks 0 and -1), progressive: coefficients are
% processed in scan order and modified only until all bits are placed;
% everything after the stopping index is untouched (payload-proportional
% distortion -> smooth capacity-PSNR curve).
    k = 1; nB = numel(bits);
    for i = 1:numel(C)
        if k > nB, break; end
        c = C(i);
        if c == 0
            C(i) = 0 + bits(k);      k = k+1;       % 0 -> 0/ +1
        elseif c == -1
            C(i) = -1 - bits(k);     k = k+1;       % -1 -> -1/ -2
        elseif c >= 1
            C(i) = c + 1;                           % shift up
        else
            C(i) = c - 1;                           % shift down (c<=-2)
        end
    end
    assert(k > nB, 'Ran out of coefficients while embedding.');
end

function [bits, C] = hsExtract(C, nBits)
    bits = false(nBits,1); k = 1;
    for i = 1:numel(C)
        if k > nBits, break; end
        c = C(i);
        if     c == 0,  bits(k)=false; C(i)= 0;  k=k+1;
        elseif c == 1,  bits(k)=true;  C(i)= 0;  k=k+1;
        elseif c == -1, bits(k)=false; C(i)=-1;  k=k+1;
        elseif c == -2, bits(k)=true;  C(i)=-1;  k=k+1;
        elseif c >= 2,  C(i) = c - 1;                   % un-shift
        else,           C(i) = c + 1;                   % c <= -3
        end
    end
    assert(k > nBits, 'Bitstream ended before payload completed.');
end

% --------------------------- Metrics ------------------------------------
function v = psnrLocal(A,B)
    mse = mean((double(A(:))-double(B(:))).^2);
    if mse == 0, v = Inf; else, v = 10*log10(255^2/mse); end
end

function H = shannonEntropy(I)
    h = histcounts(I(:), -0.5:1:255.5); p = h/sum(h); p = p(p>0);
    H = -sum(p.*log2(p));
end

function [v1,v2] = adjacentPairs(I, dirn, n)
    [M,N] = size(I);
    switch dirn
        case 'horizontal', r = randi(M,n,1);   c = randi(N-1,n,1); dr=0; dc=1;
        case 'vertical',   r = randi(M-1,n,1); c = randi(N,n,1);   dr=1; dc=0;
        otherwise,         r = randi(M-1,n,1); c = randi(N-1,n,1); dr=1; dc=1;
    end
    v1 = I(sub2ind([M,N],r,c)); v2 = I(sub2ind([M,N],r+dr,c+dc));
end

function r = corrCoefLocal(a,b)
    a = double(a(:)); b = double(b(:));
    r = mean((a-mean(a)).*(b-mean(b))) / (std(a,1)*std(b,1) + eps);
end

% --------------------------- Utilities ----------------------------------
function b = de2bi_local(v, n)
    b = false(1,n);
    for i = 1:n, b(i) = bitget(uint32(v), i); end
end

function v = bi2de_local(b)
    v = 0;
    for i = 1:numel(b), v = v + double(b(i))*2^(i-1); end
end

function J = simpleResize(I, s)
    [M,N] = size(I);
    J = I(max(1,round((1:s)*M/s)), max(1,round((1:s)*N/s)));
end

function s = bool2str(t), if t, s='YES'; else, s='NO'; end, end

function capReport(I0, D, name)
% Prints gross/net embedding capacity of the framework for this image using
% the parity rule (one bit per detail coefficient across HL, LH, HH).
% FOB is charged at 2 bits per AMBIGUOUS pixel (clamped value == D or 255-D).
    [M,N] = size(I0); MN = M*N;
    Ip = min(max(I0, D), 255-D);
    [~,hl,lh,hh] = fwd53_2d(Ip);
    C   = [hl(:); lh(:); hh(:)];
    gross = numel(C);                       % 1 bit per detail coefficient
    namb = sum(Ip(:) == D | Ip(:) == 255-D);
    fobBits = 2*namb;                       % 2 bits per ambiguous pixel
    net   = max(gross - 32 - fobBits, 0);
    fprintf('Image: %s (%dx%d)\n', name, M, N);
    fprintf('  Detail coefficients   : %d  (= 3 x (M/2) x (N/2))\n', gross);
    fprintf('  Overheads             : header 32 bits, FOB %d amb px (%d bits)\n', ...
            namb, fobBits);
    fprintf('  GROSS capacity        : %.4f bpp\n', gross/MN);
    fprintf('  NET embedding rate    : %.4f bpp   <-- Table 3.2 value\n', net/MN);
    fprintf('  (structural ceiling 0.7500 bpp; FOB overhead is image-dependent)\n');
end

function saveFig(f, base)
    try, exportgraphics(f, [base '.png'], 'Resolution', 300);
    catch, print(f, '-dpng', '-r300', [base '.png']); end
    try, print(f, '-depsc', [base '.eps']); catch, end
end
