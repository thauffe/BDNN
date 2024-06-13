#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.007724972690916, 29.93324512800868, 30.527830103482007, 31.205828853291155, 30.687847518355657, 33.23325763207813, 32.31531757337157, 32.64933793734049, 29.5085401742578, 28.468563912842612, 26.328954340455795, 28.047908567771348, 27.59920901205028, 24.083405573471566, 26.316487421097055, 24.616953649320898, 25.405223090149157, 26.526456552384758]),
array([28.49080824266506, 32.82257228339923, 32.79148034409677, 31.40255085122321, 30.52316432699543, 24.579964277015748, 25.355709386054123, 23.8348929755095, 23.98791749368514, 19.690687449424534, 19.686770978178718, 16.176928968744207, 19.30859467762352, 19.480182968849597, 17.43620319401889, 15.358179432984201, 12.992884003571183, 13.818118070730975, 11.888982974591595, 11.806213241155682, 12.885593218688697, 11.873425167838112, 12.736377790751579, 13.67095330759076, 10.391102027641178, 9.03542094775364, 11.202801162605857, 8.326142932276438, 7.182886447994438, 6.442361833534609, 5.13301283997967, 3.9111580953750424, 4.194917379869152, 1.6264746876564795, 0.0]),
array([28.557169034073244, 31.11262802955805]),
array([30.697445536596106]),
array([29.670709230877787, 28.692783499703967, 28.91552304772094, 27.533031352455687, 27.81088168751397]),
array([24.02772673992312, 27.103578950342012, 23.243731616865553, 26.99057921354608, 26.8207403081569, 26.81750361994987, 27.365472344787158, 27.419854663750336, 24.182364057634942, 17.674417565346147, 19.109321254573644, 17.297115519035053, 18.637988213616076, 16.139543134068543, 20.110220100406558, 17.476565683074316, 18.271794441986412, 14.966213309449222, 15.44115064084419, 13.534911164515607, 12.823317814321726, 12.679628939278933, 12.593427783279806, 11.677758522604233, 12.438688979025901, 12.071989834512264, 12.487854418845904, 12.888366653257007, 12.246565915552708, 12.140993574381723, 12.31226107843851, 9.166061243887668, 9.585920023282572, 8.206675733285465, 11.354283216055222, 8.029137723524661, 11.021991171340005, 9.49618632863231, 11.15699517253478, 9.296537518521447, 9.044387500800882, 9.185612615613195, 10.134715374212462, 10.466325514593285, 8.56761317912759, 11.062600272383268, 7.87601112801496, 10.707999828072193]),
array([26.47161846208863, 26.616088436588367, 24.49041001870649, 26.50014931762588, 21.505307609161207, 22.252263007464187, 20.765092724594, 16.393339828730838, 18.976510489491872, 19.07507450808961, 16.39993596510832, 16.337121882245935, 16.9873287278345, 19.344640754967585, 15.684340336936982, 14.129858481443762, 13.085048166341094, 13.453154153510013, 13.0254240552008, 13.742408198661956, 11.688220259464572, 12.92126877231606, 12.715986553466967, 13.02071655695485, 12.150099122279814, 12.91031666998541, 11.290382326963362, 11.314069335401433]),
array([20.609908024720177, 22.809461649737898, 19.758822709037098, 17.446445104636986]),
array([26.334999520757872]),
array([26.508721973985658, 25.701124437738553]),
array([23.5446725649284, 24.36810025919714, 22.679426668100074, 22.731395608653383, 22.46615221507817, 22.219029081996954, 23.018729616925647]),
array([21.56765965554183]),
array([19.680125392255636, 18.469761937785293, 19.47721133963099, 18.163182261029213]),
array([19.50888162940438, 18.973370573417387, 19.229946089715643, 19.58582980215169]),
array([18.47600258905318, 17.643737472967718, 17.732121280324833]),
array([17.05846022508673, 16.974525601136996, 17.43240161042492]),
array([16.973409809908222, 16.60399508981343]),
array([16.179562813757403, 17.302995011649987, 17.781418828875566, 16.732293782387256, 17.15649459792694]),
array([17.6088574327495, 16.920437890745035, 17.77815220959265, 16.569592927858185, 17.758830315909982, 16.67225044309162]),
array([16.626166983696653]),
array([16.17166287242943, 16.446328659691517, 13.396859153551791, 12.321477853285755, 13.681923770618994, 12.802843797351056, 13.311595731685014, 13.142133947625997, 13.065623910897697, 13.464473683349448, 13.815894219283631]),
array([16.36481934140433, 16.05728949879063]),
array([16.34806792975332, 13.91067647235204]),
array([15.988064290034831]),
array([12.921755085783746, 13.308066140718141, 11.882168627780713, 12.647052637665508, 13.348520823360683, 12.202927869760284, 12.521012513589511, 11.733634255746361, 12.676131286065036, 10.76659658977573, 10.172014303321877, 11.04427660390095, 7.929181962326851, 9.613560811308, 8.420281785414877, 8.826763513463108, 9.721469571156037, 7.070985058133619, 5.403229205274753, 6.6784908639982925, 6.4816544890683545, 4.645248573199719, 2.7241423951782924, 2.631332517456382, 2.5535968273138336]),
array([13.489980987272814, 13.576663972073312, 13.36579443323361, 11.699045969734886, 13.11967587551153, 12.00050580465132, 12.210907848712862, 11.862587549157379, 12.929742184483024, 12.386595277578126, 13.188787003505993, 13.66232136943427, 12.941000107697878, 11.706332061021394, 11.978720901586943, 12.659182739701981, 9.79590433057368, 7.279969221579146, 9.870475344518214, 11.504552675067805, 10.174805245340226, 11.492713668143123, 8.196473649935468, 10.896521735287411, 8.157555316360952, 11.39941849001423, 10.127851741630778, 7.466750216543856, 11.427005576331789, 7.796990117965294, 11.073009018212328, 8.464262087281957, 9.192101528712467, 9.195677301216307, 6.921585264473412, 6.064240865258789, 5.581584298841043, 6.375363009708206, 6.621878032068241, 5.800822575039604, 5.769949954758914, 5.454325385329774, 6.184942628271729, 6.565196787126906, 5.635842088299194, 6.508779036212094]),
array([13.894669776428788, 12.295737865277921, 11.637339590324675, 13.401004950221305, 12.41291087369559, 11.635815155798104, 11.951969592501703, 13.21215134717341, 12.128146847406716, 11.683074711865869, 11.768278894441194, 12.604531367357763, 12.139043193204133, 12.442645654907293, 13.320997142603488, 13.328799959947727, 13.724773683904111, 13.390524922101518, 13.73714990640039, 12.775583036513355, 12.603423158067038, 11.634837928166366, 13.35945348646088, 11.974233619804645, 12.059592560962566, 13.682621225541748, 10.347302158324876, 10.103166777093406, 11.308536183026925, 10.812672782354833, 10.390300631524601, 8.813899933789159, 10.634288959589638, 10.617769684935949, 11.269693920903272, 9.763074994956758, 10.089705730497387, 11.436753386195424, 9.124084784912606, 10.344876494410626, 10.143100103081137, 9.315182809842533]),
array([11.873665955097234, 12.06313187033674, 11.928923253056778, 13.766592717421714, 13.006422591701632, 13.800877800127354, 13.570249993004245]),
array([13.971838055784989, 13.609349035423865, 12.916349763713118, 13.037875261647969]),
array([14.207071025032121, 14.216427879515463, 13.744512488517538]),
array([14.179737370970948, 14.112693792456062, 12.951167761343097, 12.8439514474845, 13.235122263733066, 13.011930174874593, 13.591303459920125, 13.383619875693176, 13.043132370786271]),
array([13.45175810729486]),
array([13.330136513985904, 13.051823180641179, 13.452079777641114, 13.678497176749982]),
array([12.583063160474612, 13.064932063699603, 12.805619969253684, 12.654775836529131]),
array([13.537866701111229, 13.268867993678546, 12.874921568292084]),
array([12.587246483952597, 12.821360750254719, 12.379726435285397, 12.2303834796563, 13.461081466236804, 12.122986865410368, 12.848072617872969, 12.536408301246878, 12.368950304065214, 12.651401509795367, 12.353200498481941]),
array([13.305990074369948, 13.350768192498851, 13.389216323221895, 13.161053870671637]),
array([12.428976339761604, 12.04017731152264, 12.073308543502163, 12.290935226728653, 12.302689999682187, 12.913799433138653, 11.752682415369673, 11.797326288390794, 12.52935484401214, 12.858340341527331, 12.620455524481212, 12.42890120702769, 12.16698654174935, 12.56554493751978, 12.7902782538939, 12.166756931562833, 9.39950643551536, 9.527285316628152, 11.008555721377899, 10.196693297380419, 7.015000064680829, 6.5825901119748265, 5.519454943338312, 4.838293469883124, 4.226575153381457, 4.842914654303498, 4.147913737745506, 4.297523015944178, 5.257193391384373, 4.090489190488877, 2.792003733951965, 2.1271933343206793, 0.0]),
array([12.91630286347141, 12.980745894143338, 11.654591770010505, 12.802695357110125, 12.7252320208597, 12.88581696019458, 12.716311776731063, 11.656535622266167, 11.455205223117494, 10.573543250665631, 10.232365247305996, 10.165697167121907, 11.12895796704887, 10.55270472901098, 9.951800901249333, 10.518847088988778, 10.133070477884395, 10.65039033560666]),
array([12.008038108900523, 12.668714524907529, 12.858070932668296, 11.784615709074355, 12.623569054478883, 12.41750524463314, 11.534076114693866, 10.978703442762722]),
array([12.607431373110785, 12.325119993694624]),
array([11.930245413804641, 11.904839542616461, 11.640159403296284, 11.497506754620106]),
array([12.384853542602155]),
array([11.863028580478849, 11.013047708057362, 7.310093916100731, 9.968093483166635, 8.124404005141818, 7.51953932045384, 11.02945520843869, 10.27236617503264]),
array([11.825064019254016, 11.915633899541666, 11.667932726687654, 12.22419645636265, 11.917118011885625, 11.96108259737944, 11.935512446672902, 11.418270251454356, 8.959046860207033, 7.382222791943574, 7.833603083463439, 11.018026852473605, 8.80749260002187, 10.615238129540476, 7.9385641153139765, 8.728999148561023, 10.84426292322692, 7.939543868115557, 9.868931625235005, 10.631207616671972, 10.154060383679408, 9.730668508160615, 10.642659296785672, 9.83372081300772, 7.798397758541682, 9.464650092392215, 7.918469935584038, 9.144153345792821, 7.972621285996, 7.242255578219989, 5.889281934079362, 6.56497777888592, 6.862519850815782, 6.03420240089518, 6.782223871231484, 6.564440835697065, 6.600682522091194, 6.07503157754795, 6.46002982078715, 5.345173950861941, 6.8641299771787105, 7.226181936020599]),
array([11.903585486326572, 11.73553879239241, 12.274290489354982, 12.357054111259837, 7.939223893038104, 10.120850901553439, 9.172312197039176, 10.43761269592471, 11.302114112711678, 8.147878961981352, 7.317090908063437, 9.269828898226606, 11.444686619253627, 9.920015615102468, 9.41979907041495, 10.127243877766075, 10.193635270434289, 7.428236068508827, 11.007984889217541, 7.642146511373686, 10.173156766824732, 9.032294869544842, 8.662188006042818, 9.955684391633696, 9.449670706015626, 11.208448862859159, 7.813084193362598, 9.721963797269598, 8.465438160749862, 6.043984424344199, 5.801080678717104, 6.259602365959623, 6.2531741694027545, 5.611622479204326, 5.709273837101838, 5.654589937083893, 5.732062663415679, 6.916064675514604, 4.577782073972621, 5.1612511236474585, 4.47237926354926, 4.832305108711582, 4.508235640743008, 4.072263020012168]),
array([12.010901885244788, 11.51187222362805, 11.313122219730074, 11.248237919482861]),
array([8.45410320828796, 7.278824504272069, 5.741151177297489]),
array([11.885564361822734]),
array([11.529144204171732, 11.004939405301998]),
array([10.121919608246188, 8.823384702855087, 10.33739922778627, 8.444402800056038, 8.986319585585875, 9.930227790311314, 8.39312608666665, 11.089360702585553, 7.604937002702957, 6.122776938672215, 7.171926424069559, 6.754796976593088]),
array([11.138475657140232, 10.698320439177762, 10.62425766197401, 10.808677447059095, 10.347083174597712, 10.871237464301313, 10.585939965169723, 10.583328799241094]),
array([11.479530105273767, 11.179879495903315, 11.404501424778026, 11.051818408530252, 10.885977930009016, 10.490684159288193]),
array([11.376238974806336]),
array([9.856716958874452, 8.980294250985143, 9.127201318045941, 8.992007076670488, 9.97225351900961, 11.124847408659473, 9.01626568415387, 11.331223566086226, 11.161948254679844, 9.146351758545881, 10.445014130565948, 10.59347254280416]),
array([9.658210970790002, 9.01596329270025, 7.637557553108145, 8.278735169345858, 7.180490231090814]),
array([8.041335985743421, 9.862385269211837, 10.51102232680093, 10.32236782105975]),
array([5.320036934968457]),
array([10.603769180500002]),
array([8.537447507200174, 9.22271491812953, 9.762956040049488, 10.551393738654756, 10.623843240710428, 7.133604610708738]),
array([10.591978535724508, 10.495338844712196, 9.382667865448077, 8.481594526209953, 9.368461656585904, 7.210585942608523]),
array([9.078654466285135, 7.328921510015315, 9.682080779170395, 9.340823230398996, 8.893969629728922, 8.791249644310675, 7.401433261890359, 9.632520233449188, 7.419038902842627, 10.24503240296918, 8.252808339255209, 7.772510506149265, 7.184069235573795, 7.041990175151303, 6.67135735922027, 5.651438166055074, 6.306864457498726, 6.45922395669204, 7.213318487997597, 4.272493274398753, 5.060533842525307, 1.9167543348216292, 1.976847147533832, 1.8418024809010851, 0.6445208541301439, 0.09658252453632703, 0.0]),
array([6.9075528708781455, 4.037777620311481, 2.607014348245547, 0.0]),
array([7.762174649169941, 9.129962611118545, 8.028579593689235]),
array([8.97126216223156, 9.000544362115667, 10.17716114616617, 7.723580934571304, 7.631007533496366, 10.291718864787176, 7.209103402305617, 6.223189121175126, 6.730524762795881, 4.66239462257023, 3.0583600630951135, 1.8654078973145656]),
array([9.078708873069615, 8.933662964450857, 8.014455269309217, 9.476468579516336, 8.515202671668405, 8.967108112759007, 8.013059555003819, 9.840919942966632, 9.034663098345924, 8.98486784600509, 9.000894721752235, 9.627275841045993, 7.271367877340964, 8.437420042179609, 8.127788415656614, 7.602590952709818, 9.667532031574442, 9.768064754132016, 8.330168591597262, 9.5553375973218, 6.584455690222382, 6.489223272672938, 7.237668964104724, 6.475800317096913, 6.268563670018165, 6.06013639073997, 5.486854136955717, 6.269894972467463, 5.807987340884135, 6.723227328360647, 6.377102340018238, 5.587092431028468, 5.8027671479324745, 5.80710420379467, 5.489236512988065, 6.705639922198296, 5.012809566723793, 4.689157544323678, 4.5674741215772805, 4.821631768499643, 5.005543424271147, 4.491841145459326, 3.0223493075871763, 2.3108497698170365]),
array([8.135325128687882, 7.799784806909365, 8.653199777532004, 7.60699023531255, 7.512735508962482, 8.705810916315492, 8.786142380567977, 8.097030777368657, 9.242652537475376, 8.571659173230431, 8.135580481858383, 9.927349663255729, 7.914972324231717, 9.557199925925042, 9.710526591285262, 7.984019658762028, 7.741630460360005, 9.243519530221048, 8.716442615670273, 9.476374149748438, 8.044355707491718, 5.8078412918116955, 5.808972163882191, 5.529509559965759, 6.1649723608833025, 5.46056718466375, 7.146007405830811, 5.894890829182633, 5.504403386108931, 5.683030186558522, 6.326907906630035, 6.712091636580334, 6.368646266494318, 4.66486796193491, 4.283000274156326, 4.622380199255587, 5.1664944541473465, 4.13293628933256, 4.555699810381995, 4.7456270980286135, 3.629170271856489, 4.451371158179583, 5.312327771722536, 4.850554713892996, 4.395172736370098, 4.763630269164225, 4.517071835923073, 4.01927881690867, 3.2185219732405064, 1.8806290368046705, 1.388652915306872, 0.6806732883057641, 0.5198440087801086, 0.004145564054702441, 0.12026178173538754, 0.08100601211887488, 0.03343199082490733, 0.01039466621771487, 0.0]),
array([9.285402828232762, 9.368262071921716]),
array([7.297038587519303]),
array([9.254104064284318, 9.517766686731047, 9.568095977131728, 8.655000220546972]),
array([8.614096992977199, 8.707806390718055]),
array([8.340008899513357, 9.016983058697324, 9.069007471579377, 8.561086064082083, 9.725995054772671, 9.161867844837523, 7.874267227245003, 9.579418083229521, 8.851249847533682, 8.125217128878777, 9.534502242950884, 8.741404026201069, 8.1066263340476, 9.286256682679445, 8.158579794857497, 9.415005215571918, 7.786352697146002]),
array([9.572330349285963, 8.13971792608399, 8.137345017764885, 7.69797724305763, 9.12303893890306, 9.346982191776958, 8.069885132257888, 8.702088300441964, 8.015400621712713, 8.426815700479347, 8.342535569189932, 8.255873974784706, 9.542882461730434, 8.167323691188574, 6.816866833943907, 6.914069952140108, 6.455986751850441, 6.630550623795999, 6.853363549953761, 6.047812478355775, 3.7480603298806288, 4.984777415756848, 4.602444592741538, 5.22279968076846, 4.231714044543706, 5.117001678177165, 4.458640092096358, 5.132825921522367, 4.426960957237764, 4.511503597182147, 3.005688263680531]),
array([9.382559826722819]),
array([8.876395509449807]),
array([7.52108215783416, 7.7202442367611095, 7.962510441475171, 8.139551181535248, 7.617279929562307, 7.388874978765749, 7.399700931697611, 8.609177034434456]),
array([8.1280602810594, 8.673706340676137, 7.853406213141733]),
array([7.987788204718617, 8.265659287140736, 7.783549163886749, 6.515068170451613, 5.477564928311429, 5.457319259901929, 5.360785939154647, 6.466261267858601, 6.679755715083558, 6.252125531658654, 7.194551239085995, 6.747140214178493, 4.330233377931614, 4.292097728916454, 4.693953198497985]),
array([7.552694998410143, 7.918346418838809, 7.996225432166809, 7.331353621419944, 8.150466426005025, 6.00035812985465, 7.195688685978761, 6.878599566211134, 6.686530693561554, 6.280786414760111, 7.135393243234391, 6.265644016971557, 6.322546384870414, 5.393553082094056, 5.57858603480553, 5.996204069537073, 6.642116744324561, 5.844539833720325, 7.129269423254727, 5.690441880066086, 4.562616276300147, 3.90554286254529, 5.184765402323805, 3.855180768081116, 5.174811491562066, 5.144391915860217, 4.451357688522727, 4.539216187233368, 4.65121443955923, 4.3579218828551145, 4.387787025038134, 5.124459666668594]),
array([7.336153810237103]),
array([7.39982520056348, 7.909838895184585, 5.436869735625606, 6.738178236663526, 5.795829302650961, 6.343763533394443]),
array([6.768200781302543, 7.077910604922432, 6.46977546125558]),
array([6.0867109826619945, 7.160663764935383, 5.705964168656021, 6.942896740815247, 6.534936518975785, 7.165548917163772, 5.856211083673821, 5.939907178878001, 6.834321260489164, 6.2997431935588315]),
array([6.287729178133097]),
array([6.375132021876469, 7.055014102896705, 6.563905455728613, 6.63504903865621, 6.556321467469735, 6.173192145243885, 6.169547088268091, 6.296914392883143, 6.050111515162053, 6.371390874807959, 6.519403731506338, 6.207433458206473, 7.062242474081414, 6.950027038652759, 6.979706859470638]),
array([6.449771207175774]),
array([5.335713769286779, 6.543188063120528, 5.49056375933178, 6.749113885411573, 6.002346675229019, 5.784858253431839, 4.411913797363634, 4.5265668707801, 5.005514298811264, 5.266277028218604, 4.317386549159783]),
array([5.780429035029157, 6.3643146075633386, 6.433445707700875, 3.9633514763518196]),
array([6.441492792572283]),
array([6.233366266550182, 6.143547667747869, 6.267226456236776]),
array([5.99245542580054, 4.172803210137797, 4.177710568480408, 4.452002951617711]),
array([6.056402339988562, 6.007254657598708, 3.054900061120283, 0.0]),
array([5.8508147138162006, 5.82624679290643, 5.4951745713311535, 5.464824539087353, 4.844466334068633]),
array([5.898081444267805, 5.972374559760764]),
array([5.378818400065559, 5.096112016311717, 3.9501265290797694, 3.2164213127937153]),
array([5.544649253256062, 5.724626995680308, 4.271831744498284, 5.1998741145331495]),
array([5.040992527487462, 4.343198757100049]),
array([5.081715199596988, 3.8712829158376243, 3.722627953942008, 3.8223959854382414, 2.636205291394333]),
array([5.257342545433807, 4.81396583653142, 4.186847408949107, 4.515027612316057, 3.5596824919596104, 2.9983277123639978, 3.018911748644202, 2.2319901291814985, 1.312210014793495, 1.7634106844889073]),
array([3.6337183957834807, 4.2648386613425275, 4.634856442543204, 3.9104219010279726, 4.987459355290475, 4.132343110318214, 3.491724695453253, 2.410205203502436, 2.287380254991783, 0.9698624793683945, 1.7486460692216388, 1.7056794770402584, 0.0]),
array([4.228660300917998]),
array([4.849398623424299, 4.630044117575281]),
array([4.350127373517546, 4.88815333929045, 2.8183554719029784, 0.10295729657980777, 0.0]),
array([4.247444881733708]),
array([4.317169937949916, 4.6089050361879, 4.027644588445068, 4.704575487233235, 4.893561463413955, 4.288413998846721, 4.51267766806834, 4.643749552430259, 4.536035026321887, 1.8978099918810227, 2.3679933844839387, 1.2865785442536977, 0.8523853412321116, 1.0193225533029147, 0.7125620267636953, 0.08821658276030528, 0.0]),
array([4.739674627874711, 4.642118754154339, 3.824932361644012, 4.541868394176581, 4.668673119085659, 2.8762004020672087, 0.08065776582015177, 0.0]),
array([4.432674148545376, 4.2233739981230904, 4.70661967109146, 4.132937382554843, 4.3513349282830625, 3.990735293528233, 4.203221047755261]),
array([3.273033970650704, 1.0350554908896288, 0.5961692191033889, 0.0]),
array([4.188110792428942, 3.3664542413727276]),
array([4.468036955600407]),
array([3.7947708385436543, 0.06606861657015375, 0.0]),
array([0.10289080555501198, 0.0]),
array([3.677745015364831, 3.875166978926319, 3.2396577378856826, 1.5323014353397113, 1.7788108595946592, 0.0]),
array([1.8099703449089475, 1.2319829579162442, 0.0]),
array([2.8926898475394083, 1.628189011319832, 0.1667750524881474, 0.0]),
array([2.8733719125083024, 2.912663921142218]),
array([1.1309588244658446, 1.083940327392531, 1.0654784750095563, 0.10658652763867252, 0.0]),
array([2.8746875540599346, 1.9552291805806354, 1.5283635708082202, 1.5250258917623132]),
array([1.050763984553206, 0.0348066340045239, 0.0]),
array([0.5982617340103856, 0.0]),
array([2.075508164138605, 1.6014246654205069]),
array([2.4779274797018696, 0.8417360827917806, 0.19065655905930212, 0.12055932057839654, 0.0]),
array([2.3441007089240204, 0.8482505786494569, 1.339427522556191, 0.6157647741900893, 0.3835125072811451, 0.43786018877161575, 0.010813656405035718, 0.005192319891665104, 0.039464390401118404, 0.0]),
array([1.743902404704927, 0.5795781260465044, 0.0]),
array([1.9971088937747086, 2.0851581440995233, 1.138524703674879, 1.0383105078271204, 1.2418310110586717, 0.9646567026649587, 1.333820808429662, 1.4573969467625176, 1.6538390191052321, 0.896102924060481, 1.1797663543380839, 0.0]),
array([1.8819762731274243, 1.8172709015461586, 1.20850198912811, 1.2026937688799433, 1.2900518565356451, 1.6816445315058606, 1.698941089488727, 1.0442875163762508, 1.2320382841760094, 1.2237864971549408, 0.6871713988539125, 0.3790848788529484, 0.06544868220322564, 0.0]),
array([1.4593760148585357, 0.0]),
array([1.9105710231083717, 1.8962027164952946, 1.4175905807954878, 1.7355317696787738, 1.2054458116251328, 0.5520693611997336, 0.39044561450448356, 0.6691587789678753]),
array([1.3984216882231957, 0.9794911064989721]),
array([1.2035041539641664, 1.5841167826496005, 0.4331937424223502]),
array([1.417205182236118, 0.9747683472535864, 0.1314440870762551, 0.057349014182897354, 0.0]),
array([1.0325332216137857, 0.0]),
array([0.30507094474889046, 0.6608898032002194, 0.0]),
array([1.1578703807984931, 0.0]),
array([0.47261818818027784, 0.01995968282508434, 0.0754409066930224, 0.0]),
array([0.7172820257668138, 0.2196290815940043, 0.0]),
array([0.9803396383695133, 0.0]),
array([0.5711742955286849, 0.47041476863642157]),
array([0.5090968631700774, 0.0]),
array([0.035904860030213426, 0.0]),
array([0.5182478071162134, 0.4438499597643078, 0.0]),
array([0.6189999386820149, 0.0]),
array([0.36508346210164844, 0.01871484778617942, 0.0]),
array([0.10045178818202781, 0.0]),
array([0.1137161526409135, 0.09868290283645538, 0.0]),
array([0.09260123307524071, 0.0])
]
d = [data_1]
names = ["54"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T14', 'T17', 'T19', 'T20', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T33', 'T34', 'T36', 'T38', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T46', 'T47', 'T48', 'T49', 'T50', 'T52', 'T53', 'T54', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T68', 'T70', 'T71', 'T72', 'T73', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T89', 'T91', 'T92', 'T93', 'T95', 'T96', 'T98', 'T100', 'T101', 'T102', 'T104', 'T105', 'T107', 'T108', 'T110', 'T113', 'T115', 'T116', 'T118', 'T119', 'T120', 'T121', 'T123', 'T124', 'T125', 'T126', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'T138', 'T139', 'T142', 'T143', 'T145', 'T147', 'T148', 'T152', 'T153', 'T155', 'T158', 'T160', 'T161', 'T162', 'T166', 'T167', 'T170', 'T171', 'T172', 'T173', 'T174', 'T176', 'T178', 'T179', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T193', 'T195', 'T200', 'T206']
def get_taxa_names(): return taxa_names