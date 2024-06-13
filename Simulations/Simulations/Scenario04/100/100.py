#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.40285797573469, 32.3279238190346, 28.8355164432474, 32.04236124419425, 29.71250231595654, 28.132668726079892, 34.86156632882373, 34.253432222183264, 30.206275566694025, 28.854418051849798, 31.601913365941577, 34.515150866271085, 34.38978218598584, 34.26716370613963, 33.61131374009719, 28.94618836321858, 28.491764734806694, 34.0641106479358, 29.074783888481456, 29.775607851194664, 30.26541586025599, 29.090885607894432, 33.406943504896425, 31.416567943440842, 32.630999228282526, 30.454466787425094, 33.048747248320296, 31.532981607078586, 32.907385115085994, 31.18754697830317, 33.88822303032546, 31.07701214413385, 31.64860191006324, 34.120115149755684, 31.248524032312744, 33.61661479974466, 30.375890066001364, 31.14902570245032, 33.3965896093079, 31.44243944929222, 32.82655538341111, 29.158511412442316, 31.91310366608592, 34.309411951676545, 28.779014610780642, 29.09299108499351, 33.02432387000849, 33.21279808387039, 34.54203752244843, 32.808341465767604, 34.64228333282114, 26.205387130180764, 25.232574615531174, 28.05500285270348, 27.274374093811137, 23.546472093547393, 26.297953074793426, 26.52954967810091, 25.14846564120606, 25.452625366823444, 23.625949249386757, 26.71921356326116, 26.57991231211077, 26.536479396550327, 27.587908342449538, 26.6042399882048, 24.40445270434896, 26.366674013933714, 26.64963445152173, 23.251317637036514, 23.03264070375927, 26.389625758536038, 25.691133707024324, 23.893676223847592, 27.596015908351227, 25.523858766583892, 23.87244860857767, 24.195476423393593, 23.962878226742276, 27.83864161957367, 24.64568839548646, 23.21610128166447, 24.546196395178704, 26.02207838790102, 27.05030421460609, 23.252965078703085, 26.04152333176658, 21.768336689844766, 22.71341963752062, 22.537858343379746, 21.476672026226737, 22.257978941785748, 22.55685007745912, 21.727158282113745, 20.56776563975987, 22.427617561827045, 21.11577750533934, 21.057271296037708, 21.045913717623954, 20.635555873847565, 20.398887574622073, 20.06892414925528, 20.367975354763086, 20.019593973181397, 20.234533929836008, 20.233275603954485, 20.402374588424365]),
array([28.163625317671634]),
array([28.419544448691163, 28.604365103817557, 28.823488410035765, 29.05012399977471, 28.740757413483777, 28.914320740091824]),
array([24.4116249611855, 24.629022726812927, 19.576348163929932]),
array([23.683552395427583, 22.700735099389174]),
array([24.571778729028118]),
array([25.129587617464203, 26.5634425654007, 24.745925820984304, 23.11008296158636, 25.271697781940446, 23.750583631527405, 25.25849159697246, 22.418951404604403, 22.875362009381345]),
array([23.99471009742808, 25.94081163954124, 25.333223213684548]),
array([23.62728846960376, 23.835162018151482, 23.54250324857656, 23.127182626214108, 21.23477792220258, 21.792011836542002, 21.911939188799337, 21.394163231569728, 21.819563038111877, 21.05411157960578, 21.603625317345376, 20.881624258089094]),
array([23.80375173931081, 22.2686036418755, 21.455308188278703, 22.31935378605059, 21.516778393133404, 22.218252987946748]),
array([22.915958775775568]),
array([23.070676445499146, 21.880277023335278, 21.19543604910581, 22.36328567626312, 21.428439757789164, 22.714990528646272, 21.37933004200875, 22.879024375481393, 22.906653577820993, 21.704063591661818, 22.2547741736085, 22.132730409467925, 21.34852684437643, 20.313990986212417, 20.381621108889746, 19.999573459409785, 20.199870336357346, 20.079421023730376, 20.127213582907192, 20.424876611083214, 19.974695285515306]),
array([20.57237938041869, 21.692933532929302, 21.181922488431024, 21.511842365636973, 22.00711617922933, 21.289995199498644, 22.415728207443287, 20.83313831125123, 22.272494497442775, 22.717500115531422, 21.91323293599728, 22.325645457708493, 22.32104894494769, 21.366695345354685, 21.864727069190153, 20.299146242802095, 19.921221774706368, 19.937196416370636]),
array([20.677024438268887, 17.19625521850196, 17.196992558134987, 20.212813185051708, 17.784524468965984, 17.23444384596629, 20.25282790968488, 15.04619518364389, 10.85474592355538]),
array([21.096158003711306, 22.02108145767312]),
array([22.1696059837688]),
array([21.665019019462374, 21.71969066796233, 20.948810552650897, 21.379120584763836, 21.735709569581196, 21.421569383802034, 21.891377288717067, 21.530703548382593, 20.617940201561687, 19.4683033331272, 18.46635148211222, 20.344347330198516, 19.2822960517255, 20.29688120416113, 19.879716511403384, 19.03635634798911, 20.115356120803902, 19.45542702559691, 18.530908577405032, 19.261289675269516, 19.22674635957852, 18.704026911922988, 18.56690412631083, 20.15167908013676, 19.885985619598436, 20.189356125489425, 18.670558423684973, 18.78652445907106, 18.44928223781444]),
array([20.517355633949702, 16.802002817414184]),
array([20.840759213790395, 21.381114363901805, 20.99495584993299, 21.246353348211706, 20.960242378632287, 20.983430125379837, 21.065817873188884, 20.052776473061083, 17.41250545357857, 20.417293248623903, 18.690212189279485, 17.89686998785833, 17.551174931884287, 19.44752921039255, 18.34479811220936, 20.034012123391207, 18.462101023821038, 19.82413007194725, 19.013036087015678, 17.28615904018053, 20.218698545487634, 19.77318819765125, 18.7903911919417, 19.983374617741145, 17.405350343888394, 18.46872734695307, 18.940022530201176, 17.835117944543967, 18.16920123575153, 18.020188635513268, 19.383217719609604, 19.86003494173239, 20.39335761646693, 20.09208838633753, 19.552581971974124, 18.43630422982116, 17.360367522766712, 17.493964801260496, 17.30344440830583, 19.342731510666106, 18.30668996507515, 19.492768711772992, 19.949675023524325, 20.31997079374196, 18.496617836838393, 20.34987009804125, 19.197728891640963, 17.341934535878412, 19.7178547865099, 19.448340109434334, 17.777210643766626, 18.55997309829955, 20.35024415036421]),
array([20.66284998944002, 20.548530907180915, 20.60338300631858, 19.266103040746383, 20.056069245713392, 20.21359035922646, 20.419742359200054, 19.795405006708407]),
array([20.206260526813587, 19.301881134439324, 20.436032676405514]),
array([20.164119417757696]),
array([7.886211969920595, 0.0]),
array([18.449505547448, 20.24502087276164]),
array([17.713904084528167, 18.149467481336988, 19.170862017272498, 19.08929024820698, 19.87300460532256, 16.54689970473735, 19.49623496071385, 17.158082492869152, 18.557122452134905, 19.61157336042385, 19.881008803198384, 19.363231554167147, 19.917309025790903, 16.760800051593407, 17.000098633906447, 18.325139166871644, 19.6186971152946, 18.832583646292097, 18.389269869356266, 19.426246261210228, 18.167469063687495, 17.250822078149042, 17.332731343315977, 16.89933386531472]),
array([17.90610758548458, 17.69394373810174, 19.097046155164108, 17.98994321640414, 18.08965113141375, 18.897886168988002, 18.229442014654573, 18.802665161224454, 18.505371374211578]),
array([18.86449663054637, 17.55795063459903]),
array([19.16399339089208, 16.87389744036978, 12.286007267053058, 10.519604509453332, 10.386014130053843, 10.793652035667778, 9.219233601319079, 9.341716801309985, 10.977051701956363, 8.611596249280902, 6.775869182857767]),
array([18.84753845577566, 19.276707051054526, 18.897365198566078]),
array([16.715306814527832, 17.668264560763284, 18.089116498150954, 16.2527204107793, 17.67870702466366, 18.043799218118263, 18.294474514372673, 16.77927356035635, 18.49038387524058, 14.61337984552016, 12.691294546256378, 9.007816666481741, 10.243687677817135, 9.9733017952798, 10.020338141982334, 9.117053822797958, 11.590501612443912]),
array([18.01624036596755, 17.975777991261207, 17.932856237326636, 17.515311101468377, 17.70841649356161, 17.760532376549016]),
array([16.275429540468526, 16.05584541004841, 17.627351237265707, 16.65551810706598, 16.908140736744805, 16.4484611252336, 16.392987650630687, 16.163339275409015, 16.436597772142918, 17.764490150249348, 16.95927698442711, 17.82489694086038, 16.435238554129494, 17.314186715915774, 17.794060082114367, 18.14242223638412, 18.116768248635015, 17.986273573152676, 17.230515688255135, 17.962462930753297, 16.553902622919345, 17.83585480614448, 17.820801957803727, 15.764193869661785, 13.568421852861336, 13.290042111980116, 11.621846728252967, 11.204585227227145, 11.33977789187175]),
array([16.84548928004321, 15.981683990467292, 17.090850767191082, 16.625121414600272, 17.522569605925472, 15.309471505192256]),
array([18.21904095002393, 17.92566171863412, 16.03086090786493, 17.174790894073325, 16.004902303013093, 17.587353126654154, 16.29332698140748, 17.714669039979114, 17.212838951853588, 17.04850184746671, 16.42640118729425, 16.327521160588574, 17.201207344987054, 16.0689205153566, 18.014202342043536, 18.31263086072195, 17.412253009453057, 17.566950995695215, 17.7297615092599, 17.00839180987695]),
array([16.541887664005564, 16.135531188684485, 16.205357348302105, 16.408972215147738, 16.057803159762837, 16.88906738040938, 16.7135975674586, 17.444654751535882, 17.533000914834275, 15.997872421659448, 16.063666446379553, 16.093046227434527, 15.992561169432193, 16.299615763450028, 17.50262083496866, 16.764520088376248, 16.520247082599337, 16.86433064955038, 17.69205032007406, 16.216887024951262, 17.295263374492908, 16.642464175281415, 16.51954328286368, 17.177654885857844, 16.666683151423683, 17.023583828531656, 16.739057232561027, 17.521249278939667, 17.648686747323744, 17.12293981701005, 16.739126654018918, 16.81447617408767, 16.67894211532091, 16.583349509354814, 15.944920790306321]),
array([16.69730725083575, 16.907962743131097, 16.85738578967122, 16.254588156846285, 16.66897680181285, 16.46532422566169, 17.29210732657626, 16.73506488913483, 16.11352182703621, 16.773303875674525, 17.210502758830454]),
array([16.672145353065382, 12.799131037773485, 9.100166401942342, 3.9905403623661946, 3.5823704369009843, 3.1744362127315022, 0.8742101731304588, 0.6273024611271332, 0.0]),
array([16.28527806073683]),
array([12.930466011399437, 13.355290680502678, 12.991768819006948, 11.50424843795547, 7.489821238123721, 9.139298127894984, 7.743252527752643, 8.325774717643473, 8.833019158614686, 7.310199365187969, 10.792950690317118, 11.07148268105151, 8.655517464036135, 8.866120245922748, 9.17693681891197, 8.411210104212616, 9.055590872490288, 10.073271073546145, 7.343468637139869, 10.0265202980763, 7.604613795629016, 9.23886737714754, 10.088213170242808, 5.679452074394618, 6.02773033835907, 5.8727687334143575, 5.692916437455089, 5.722965913932977, 4.653975122837585, 5.288261761224151]),
array([10.141052010048206, 11.40798990583152, 10.81170111719836, 10.979633521463171]),
array([12.359930215674142, 13.418709486742037, 8.907198961147829]),
array([14.386646696670716, 14.38685593464309, 13.132713908065394, 12.855609879971816]),
array([12.135494367074719, 11.401846883898799, 10.263937398133647, 8.478603841197632, 9.784530441064064, 11.393585026995362]),
array([5.540055110620016, 3.5556901554300966, 1.1445841356695559, 0.0]),
array([13.51275360121523, 11.439610754064061, 10.210697832122175, 10.514713169123949, 10.426443192046696, 10.717562933367198, 9.949916770265682]),
array([10.603813696422128]),
array([13.330438342423802, 10.792022686708457, 11.353196980134259, 10.500897716575803, 10.75506839557275, 11.503098333875844, 11.435396631574628, 11.024826647324176]),
array([10.939336733928828, 8.226308445483607, 6.798302778136189, 6.78375629421999, 4.006765907599914, 2.099795768168673, 1.896277300500251, 0.0]),
array([11.514566410512598, 10.42066827126648, 6.57038810920952, 6.375292824480366, 5.882731824830428, 3.073908294806814, 2.9096912903451555, 0.0]),
array([10.02686607518698, 8.314178423673987, 7.00043229914483, 6.807543061066011, 6.126145214195386]),
array([12.48580661875613]),
array([12.556649251706828, 11.9089072072644, 12.044368637277362, 11.968708882202758, 11.849464733122863, 9.664647480135711, 11.000982601047289, 10.127287259677379, 8.057873976647825, 9.823416164799308, 9.500547795877933, 7.720876952552149, 8.828287700640319, 11.452810863074895, 8.705814521421193, 10.208861113656134, 11.269258905077649, 9.239865774040858, 10.659749297521895, 10.029747818236224, 10.322516188250999, 7.976301775953429, 11.28442362209907, 9.078644424010946, 9.164594808567543, 10.52596371814022, 10.94803721243419, 9.050917767258333, 8.969598218731143, 8.430250769831131, 10.28794885981372, 7.101321692916749, 5.867779584366204, 6.289846701648031, 6.779249919896995, 6.450407242540273, 6.806770066260041, 6.3236169448961475, 7.130161824405067, 5.961163583780531, 6.891851665694079, 7.216608703453549, 6.274344259909092, 6.237149707092637, 6.345031298866781, 6.57087190047132, 7.13276947438759]),
array([12.158693152237554, 12.545219208147559]),
array([8.293434036489153, 6.477492346029971, 3.30745439555631]),
array([12.191720250363188, 11.14343649124451]),
array([9.888618578789748, 7.919169062477175, 8.782314828101677, 9.161743822845118, 8.143932845831515, 10.307664622467435, 10.420352156489805, 9.196260037283125, 11.349120351088075, 9.453289211090482, 9.976497382636545, 8.309108464840168, 10.318635329799706, 11.428149098516757, 10.051633927193263, 9.558711287034193, 9.153530780174737, 7.673658646288606, 9.914661997510137, 8.141527886466946, 7.92644624799042, 8.36578233286053, 6.942498357337781, 7.1560639343231625]),
array([11.503637577942822, 11.530074386539471, 11.350066942615577]),
array([8.437140693880796, 9.901109395653624, 9.223406060651332, 9.568092600342466, 9.950527505972016, 9.409207073467716, 9.127722519967216, 8.637159072904357, 9.696830554265398, 10.726731190941308, 9.28169177160105, 10.278682871083868, 9.593808126246334, 8.639770201444136, 9.093975153416668, 7.181659545157158]),
array([8.724759093317076, 7.600457854981084, 10.616198075245189, 8.713345042155222, 8.003387903081023, 11.386225961792139, 10.234127447746932, 7.603974439857561, 9.763250975012657, 7.827846509166506, 11.316591520001495, 11.145146072101065, 10.558195234778646, 10.5682537615953, 6.357164600106831, 6.700394769398214, 6.661060792632292, 6.5514257879317706, 6.2897705486812585, 6.291966320742806, 6.104606593440653, 6.820629462255071, 6.795067253359339, 6.123357863797453, 6.936929109401276, 6.349469748339407, 7.026519160697939, 6.448243248300887, 5.440394971638778, 4.895938751688804]),
array([11.295113592299407, 11.277977870598855, 11.373990348769695]),
array([10.960023985080749]),
array([10.057528001523073, 7.9889526942635545, 10.642925988976982, 10.074107324129486, 7.454959800690636, 10.861112108148047, 9.821977579252243, 8.396303481715133, 9.09346527787559, 11.294793582536043, 10.087625191595524, 10.963668850250937, 11.28939064220382, 11.08226727558571, 9.392600074446467, 8.745649377489142, 9.424473586366815, 9.654519360805766, 7.28089951740247, 10.420395097634682, 8.916545247414835, 8.28423974974903, 10.417249613969531, 10.022463313556703, 9.411947816658293, 8.10761962177725, 9.738679624026629, 7.600214615046468, 11.225603558300437, 9.790157262912892, 8.714153076958834, 8.221528213097423, 7.956339663866074, 11.14594618731716, 8.121634370934913, 7.197924294807514, 5.453163754008184, 5.687896525108917, 6.988114545019228, 6.462146944095615, 5.893527380160585, 6.152340015518223, 5.569668685667976, 6.72747466749348, 7.060108553818289, 5.99630625970111, 7.236030910263761, 6.9569419105214845, 6.4342369700274755, 6.39180810908204, 7.091836730354608, 7.147421414128321, 6.335251587304608, 6.232143272045647, 6.302002519943701, 6.772045143079597, 6.385759136676204, 6.608176480789614, 6.938537396880995, 6.582608471106881, 4.727610054553467, 4.956232607527474]),
array([8.214571138438865, 10.189170740262421, 7.818344080894052, 8.041867023365821, 10.157482195240043, 9.86964915352879, 10.718943051232852]),
array([9.598436398990364, 8.412293862099519, 7.887766357081073, 9.627383979721554, 9.840084069544364, 8.712086549996938, 9.157546054258521, 10.572882191446494, 9.431689921762414, 10.82361575151983, 9.623958589363161, 9.828494002362257, 10.763938219242204, 8.821934570369848, 8.3750982850157, 10.35314747789548, 10.744289582358594, 10.354611795489783, 9.81064833700277, 9.363077666603596, 10.54590964939703, 10.835980451192817]),
array([10.100057720642653, 8.734512208436321, 8.879720527620915, 8.25831919817355]),
array([8.628670613134194, 10.45118719335202, 8.341304810014359, 9.314788414144536, 9.035141600362454, 8.77132469153591, 8.96457831025762, 10.16766930777172, 9.878875703122036, 9.062349458127219, 8.808316853059189, 10.027122425335874, 9.053646616482638, 9.016717380226492, 9.693178957877798, 8.28453228148485, 7.936800941934573]),
array([9.887681409959068]),
array([9.841981095090789, 9.704889491582982, 9.072251438726639, 8.501921375999508, 8.780761789127956, 8.288284554631794, 7.547221734475949, 9.414961622598318, 8.31971521583205, 10.266586722334234, 8.415288066217073, 7.5028683995326375, 7.405877986430852, 8.382467401014152, 6.8782017400502315, 6.898093555509918, 6.513453390544719, 6.520505734638121, 7.118158612193018, 5.9370907199717795, 5.824868609572958, 6.228521133227924]),
array([8.427060203120567, 8.452628275603088, 9.642611870682058, 7.469205630400843, 8.88760921613601, 9.74518620118361, 10.151464782314827, 7.020445577063983, 6.444906236286072, 6.8283503817425935, 6.850805932857792, 6.488301455995997, 6.971094149211117]),
array([9.909144105435418]),
array([7.6665785424792485, 9.089225485019636, 7.258173568957989, 6.389891854133241, 6.698210809206306, 6.549454589223732]),
array([9.573760597376289, 8.423641128841465, 8.420186323262909, 8.472453214359176]),
array([7.814587660219922, 9.308243154851676, 8.80153245364644, 8.782761340666843, 8.149809022894495, 7.9568623272186, 6.407744154324007, 6.587631514851826, 6.758044836858403, 7.221905327796176, 5.439209283590522, 6.9613408268242445, 5.8773359530848115, 5.931909900199548, 6.003212223483295]),
array([8.893049202970708, 9.481256262264498, 9.462870537013568, 9.001365078415212, 9.237282337191814]),
array([9.149137271194203, 8.092944675170667, 8.557137454923362, 6.6942010860179835, 6.1133442061254595, 6.259731620996957, 6.382976681373131]),
array([9.10915174438312, 0.0]),
array([7.655512759397556, 5.678475764931765, 6.979143228099055]),
array([9.159967744641708, 8.626204108785547]),
array([8.410855769202545, 7.939158650315745, 8.980218138295522, 8.465032383794373, 7.87143764230345, 9.133398682637257, 8.146250347620203, 7.455522274393079, 8.962877607358127, 7.581630480149713, 6.448663969608227]),
array([8.755015791725015, 6.382616190395113, 7.175822142897587, 2.346596286616815]),
array([8.492183129685692, 7.8991391322692905, 7.948054440485528, 8.846675168579798, 7.908456625752866, 9.044217285049353, 7.613853856891098, 7.852736273773302, 8.579201779576003, 8.300138819733583, 7.283145089665074, 9.10600661014225, 5.458827702517195, 6.661837662391012, 5.4089274996468975, 5.420715722137117, 5.78778915055533, 5.4539207950163195, 5.356085941068313, 6.978624436017844, 5.546446500271986, 7.0195466759508145, 6.717107742170959, 7.232667760803074, 6.92911367072948, 6.381400982066702, 5.619262430788099, 6.006013332939853]),
array([8.940851481799346, 8.509946891097874, 7.615412554364698, 8.518871117290503, 7.167906897484307, 5.922140036837238, 5.820173983280899, 6.585332746027438, 6.638234858195564, 6.722410693095311, 3.2826089369288884, 1.8186031033677406, 2.0944532502135864, 0.0]),
array([9.120703713719447, 7.447974204823275, 8.877833598842223, 6.510955860867767, 6.588918573360221, 6.512486716825271, 6.240198136736195, 5.802690379488005, 6.839527112904983, 3.875338805951081, 2.939425787855531, 2.808220116097449, 2.444163559550133, 1.0609512375503107, 0.4675232040026353]),
array([8.99507807686558, 8.62562542683232, 8.979884831277975, 8.176627313107652, 8.199313645381455, 7.828399504547161, 7.732366074365459, 7.887668291449421, 8.5994047710221, 8.860907085032899, 8.004493770236904, 8.777018480445054, 8.730028005472231, 7.746232549490065, 8.88194680718189, 8.927566656738676, 7.124520790276066, 6.862227672147764, 5.788792954517904, 6.865550806272594, 6.072974234720697, 6.443595323029367, 5.531182852633575, 7.213627239441078, 7.206489076799548, 5.5019617967929415, 5.722017599942844, 5.985756272377065, 6.64022991016098, 5.673067577537897, 6.9064175666653975, 3.949826836392447, 4.121147061349124, 5.312329823406336, 2.5842051840968043, 2.9139389832451275, 3.3218003357179646, 3.4898863065329233, 3.1356335215519775, 2.973229718143032, 2.9225331399821988, 2.8641159132711147, 2.850230473526283, 2.648067574296409, 2.2382380979457643, 2.200581247017902, 2.2681670252911936, 2.3927173684745875, 1.8951181712388838, 2.1635498267159816, 1.0251691061400872, 0.9547497864208992, 1.7216171450669184, 0.0]),
array([8.737779625366848]),
array([8.56820273340788, 9.016357262959739, 7.306954213350427, 7.2651693474133925, 7.739188282816494, 8.055541868914418, 6.9659908246966875]),
array([8.63340809813355, 7.556064934939908, 7.746216982186187, 7.392814627543785, 7.585511674639822, 6.698803917165899, 6.576852651699456, 6.04128465137751, 7.1725106219882235, 7.205243797214254, 5.992643565244215, 6.725562950950337]),
array([8.564943475906002, 7.656917824033603]),
array([7.918150147596306, 8.575425229758746, 8.465659545986613, 6.306554874989244, 5.46470244110168, 5.974801061498388, 6.905157651578314, 5.429721450658627, 4.395780488822727, 3.4129789124286822, 1.858288138030665, 2.4425546026239693]),
array([7.461585936777114, 7.685241412946298, 8.55013673468125]),
array([8.019126166697273, 7.273934557046784, 7.779887561645932, 5.68144100206151, 5.726562041262826, 7.204129568322694, 6.6293242028505635, 7.2389111513304565, 6.9686562603371245, 5.573108474942121, 6.32659648883896, 7.213378247409806, 5.732686010912488, 6.593899942991869, 7.230767045219247, 5.816296966325195, 5.637622339921654, 3.2709474700705585, 3.5957425997723087, 2.9991948951356866, 2.932927208293309, 3.5952427930229978, 2.0061941339388087, 2.52996046047034, 1.8637123621196023, 2.33998574727263, 1.1042394222585283, 1.3000199828603074, 0.5785886548456272, 0.4059432972117585, 0.7050516835167276, 0.18282402800979702, 0.0]),
array([6.061268927531531, 7.101468167834272, 6.849338754821382, 6.579191773487263, 6.0710716521548305]),
array([8.116519038048775, 8.309204884581922, 7.710413883229042, 7.275952743279163, 7.765513320429972, 6.13463357097135, 7.129141970281229, 5.927144718583177, 7.0473793167246415, 6.714444346188218, 7.137861928108546, 5.39159075015973, 6.0490244543147975, 6.329973565580006, 6.965374808009561, 7.0289230611805005, 6.089623127786345, 4.875457139358238, 3.2495436001188485, 2.9852401801329504, 2.914271041824988, 2.3189040562838095, 2.547587652517961, 2.4251221261237794, 2.3157247955747264, 1.9119969922541196, 1.8949281641077615, 1.8292507202212436, 2.459961539228115, 1.2148359930570596, 1.142955634822465, 0.0]),
array([7.720380955010971, 6.918866236006141, 6.861583000715784, 6.980507852315023, 4.814503298127318, 5.2272306674816456, 2.729369301685585, 2.2729592911760728, 2.542265202509345, 1.8671867574955336, 2.339749765507661, 0.8495106739041008, 0.41687413733120005, 0.0]),
array([7.760795666351452, 7.369846280464862, 5.444499191879216, 5.899494707897496, 6.224539457783635, 6.182890097415911, 6.072747915289642, 5.794285866557361, 6.423828763447559, 6.8749385372452645, 6.70878124106256, 3.1796070377735535]),
array([6.792004653275886, 5.483343614807253]),
array([5.164919204773604, 0.0]),
array([7.441855346668863, 7.273856014263718, 7.320325059871531, 5.664686211265514, 6.681783837021676, 6.705641484426476, 6.32331758110859, 5.489687103744667, 5.418166839578374, 6.638707860276966, 6.808970290773153, 6.209043270033112, 5.567371490121213, 7.200780655365597, 6.497901685173871, 6.096379937487114, 5.909991178310408, 6.083690138059165, 6.42966488263142, 6.659828930691505, 5.414929338136796, 7.175502310902418, 6.781206962329903, 6.0648793973563135, 6.2466167108187935, 6.032039924276806, 6.567146225217793, 6.5289172349437745, 6.808083795286938, 3.982911191574911, 4.068327291936994, 5.179342991065546, 2.77756567592977, 2.972054389970523, 3.4362773051088817, 2.9497240902500925, 3.196489040990478, 2.943230358534595, 1.9961699880305388, 2.0499100767573677, 2.455094680799162, 1.9224159117445743, 2.076047749651153, 2.1228644353719064, 2.553210803217449, 2.3985462023211594, 2.337873472529691, 2.0011265063328927, 1.9061625641411155, 2.0358163918859358, 2.1844362047320476, 1.836283455483908, 1.9995425923655483, 1.7524328616643492, 0.7314040491752377, 0.0]),
array([7.063769848165501]),
array([4.724756672216795, 4.697248490471062, 2.2832272377006024, 0.0]),
array([7.384098327434587, 6.99972266310925, 6.401880541583346, 5.934178444303889]),
array([6.8087228146770755, 6.88473260470564, 6.529013266670331, 6.6077929669249515, 6.614224082067549, 6.952342566172278, 6.902944967803023, 7.233248601640161, 7.04898109903217]),
array([6.3437788802861, 6.5346523960633665, 5.766742415642722, 2.6926124165985357, 2.260984180651829]),
array([7.102686167185438, 6.74271704122044, 6.985062039636511, 6.720025249839348]),
array([6.02439428316294, 5.690999412469618, 5.99288069794173, 5.3784449974368265, 5.738002294940414, 6.0633787554792296, 6.385996489133616, 5.554695750103039, 6.982891898245288, 6.871544363343019, 4.284928299117416, 5.200630211168196, 4.35656895368103, 3.1990302582341843, 3.065683877789515, 3.557213779972757, 2.96880336020169, 2.892551776057384, 2.9823892900524074, 2.1264621014265357, 2.145719667625027, 1.8188713964769436, 1.9143760680932727, 2.475098224914094, 2.292035711554275, 2.543455210208843, 0.6068927852288557, 0.0]),
array([5.7160052280647236, 4.094874702476908, 2.6376542346319196, 2.253732663312729, 0.0]),
array([5.86001440136631, 5.784850988043261, 2.9672126123781126, 3.2728218824361357, 2.475566969518408, 2.4279627898217067, 2.3082381033731405, 1.558707437084931, 0.0]),
array([5.751366440116925, 5.702236934811824, 4.774916657439391, 0.0]),
array([4.147972115803306, 4.977574070737664, 3.0116553958593166, 2.4318030181296098, 2.5573627743138663, 0.0]),
array([5.149424477903272]),
array([1.1686851777238865, 0.0]),
array([4.433442854666457, 4.894618134687171, 3.333295971461713, 3.232648456799605, 2.1489620887476977, 1.851133784249015, 1.876990090407065, 2.198735057193981, 2.222826394414393, 0.8344354289966568, 0.0]),
array([5.088268024191412, 4.014526714966285, 3.1129566549887873, 2.445023566527847, 2.2354443192610773, 2.3573530012177812]),
array([4.099281332600451, 4.391651107844861, 4.274613361008576, 3.3748742767057576, 3.224145765739122, 3.023437589349898, 2.6402470232694313, 1.8341448198307335, 1.8733132883069692, 2.5058075350621367, 2.0581831326929607, 1.9546918023048994, 2.08549492093015, 1.8911615024300725, 2.099852853041036, 2.5454119809805005, 1.1714530135950123, 0.2837293925641268, 0.0]),
array([4.746696450041357, 3.3568561908488603, 3.1767271535701465]),
array([4.205840767103107, 4.3620205893709505, 2.9760917344272624, 3.5156866346308315, 2.9794363657520275, 2.811543089597383, 3.5121383698198763, 3.0727072651845644, 1.9982238222630557, 1.8835233353687397, 2.326328861288437, 2.0506264382972708, 1.8890142817670237, 1.884876950254787, 2.1484347844709197, 2.3277785571197196, 2.4225958556324825, 1.6960526483049785, 0.8023432060538976, 0.9123744713602997, 1.1930146051481865, 0.0]),
array([3.1012557816929327, 1.3457943658682945, 0.0]),
array([4.055680512464245, 4.228895892463758, 4.2510272522204655, 2.1367310513497864, 2.2279183371258577, 2.5680900070344115, 2.5659202204545966, 1.8010948900841766, 2.2255593929423063, 0.0]),
array([2.9356168830523317, 3.0152042490744844, 1.8778876052808993, 0.0]),
array([3.8763831463369183, 1.1837600082184259, 0.0]),
array([3.6955687596711173, 4.086166110416344, 3.4243630777614453, 3.3622575046722876, 3.3037597321470145, 2.918748495444214, 2.8120994002738424, 2.805654701523017, 2.956378416469601, 1.9444465774690118, 2.496525544582897, 1.8608182680324568, 2.253398245562315, 2.141726092882938, 2.3475230076062856, 2.4974125464653656, 2.133827763338269, 1.9320103826889077, 2.194793091785482, 2.047671717700357, 2.163455586621849, 2.256059717556056, 2.382711885930837, 2.357834711668215, 2.4513959469009876, 1.009107807162135, 1.7900847513343197, 1.1443856323659323, 0.8948697635801325, 1.119047915984294, 1.621144860036715, 0.20152420456882147, 0.06642223336793499, 0.0]),
array([3.4282202305944716, 3.0546691090345903, 3.1171009969111143, 2.58683016580524, 2.4348639293106333, 2.403225303518048, 2.552578600179982, 2.371340657618609, 2.417956797419419, 2.072059873949609, 2.5052806618485635, 1.9240460348698571, 1.4881450754740413, 1.3412850304845567, 1.3330452562159985, 1.1600139194568273, 1.2392876982720864, 0.2626291809740182, 0.07485414409386967, 0.0]),
array([3.3137227835441334, 0.0]),
array([3.2552970784566875, 3.503647287816048, 2.1185564738802367, 0.0]),
array([3.744758433548841, 2.6340164185824317, 2.7252847052572173, 2.716269215928615, 3.1561696129253796, 2.5835847729514203, 3.1975335183989015, 2.587198544649622, 2.509104659032576, 2.184982902300204, 2.037015314293215, 2.53551118516372, 2.26455915946662, 2.543548723343183, 0.0]),
array([3.653862152788496, 3.617900776079646, 2.973434840011989, 2.868434577172634, 3.3305031452598146, 3.322617748965945, 3.1243920444726903, 2.6128461331278574, 1.9061326678477744, 2.2780539099767125, 2.3317732147433845, 1.8559516646226264, 2.36387874729134, 1.971242811374535, 2.020064013954868, 1.9827545063006005, 1.9578234908986039, 2.186287870688758, 2.3855709793368605, 1.101803756492492, 1.3263061147563384, 1.6612811061016697, 1.7367418972123811, 1.103353855485012, 1.3837697075308582, 0.9215312237521396, 0.0]),
array([3.115400808822794, 2.666435334959181, 2.7083289844969913, 2.135666746941441, 2.0046921544043173, 2.011315642873854, 2.5366867113604465, 2.141715261325422, 1.7889890495471494, 1.4511205885200256, 0.30281574595882615, 0.7361576263079117, 0.0]),
array([3.3665368445986883, 0.1915672650503667, 0.17686249306425217, 0.0]),
array([2.1091860390130295, 0.0]),
array([2.4232948757596886, 2.0417229317832524, 0.0]),
array([2.714286855361463, 3.2621613405760472, 1.9898053173496377, 1.9788459955435713, 0.5403795810249359, 0.0]),
array([2.6378336745972364, 1.9503292251672788, 2.0370694320796146, 2.5707871179170865]),
array([2.590802666808112, 2.117471120602771, 1.9968748687036038, 1.755647349157158, 0.0]),
array([2.9619607905526912, 2.6085528384353176, 2.408360403262508, 2.56367526681863, 0.8795070145581722, 0.9098033402182484, 0.0]),
array([2.3197033225047012, 0.0]),
array([1.8758604048750587, 2.08727615530695, 2.0706686523889406, 2.2697526841972184, 2.0465032053340417, 2.385118998743708, 2.201959025517539, 1.8316041120545208, 1.569115979720722, 1.612560938223464]),
array([1.8711859462632128, 2.236865217242654, 0.8310282150687343, 1.7492081999274431]),
array([2.106382274985997]),
array([2.0413496361755104, 0.8561245495044145, 0.0]),
array([1.87282596878674, 0.7960056529443296, 1.3973739699123622, 0.7852660660115991, 0.138626111416276, 0.5666503879182799, 0.2018660951272938, 0.0]),
array([1.5887243155642923, 0.0]),
array([0.9428835085234881, 1.0989575046602713, 0.0]),
array([0.22950463794181164, 0.0]),
array([0.36784338591962384, 0.0]),
array([0.30605244589998226, 0.0]),
array([0.20110314700172058, 0.0]),
array([0.1263712575079031, 0.0])
]
d = [data_1]
names = ["100"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T2', 'T3', 'T4', 'T5', 'T7', 'T8', 'T10', 'T11', 'T13', 'T14', 'T15', 'T16', 'T18', 'T19', 'T20', 'T22', 'T24', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T44', 'T46', 'T47', 'T48', 'T50', 'T51', 'T56', 'T57', 'T59', 'T60', 'T61', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T74', 'T76', 'T77', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T87', 'T89', 'T90', 'T91', 'T92', 'T93', 'T95', 'T96', 'T97', 'T98', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T118', 'T119', 'T120', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T131', 'T132', 'T133', 'T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T147', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T157', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T166', 'T167', 'T171', 'T172', 'T173', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T187', 'T199', 'T200', 'T201', 'T202', 'T204']
def get_taxa_names(): return taxa_names