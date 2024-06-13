#!/usr/bin/env python
from numpy import *
data_1 = [
array([25.796002855240342, 27.97118220943716, 22.749600472519766]),
array([32.47339883300753, 33.29768711508742]),
array([30.20130623825449, 26.343857549692967]),
array([28.243241639478438, 28.19021255978853, 28.145598693755296, 25.612739095734998, 23.845770150810008, 27.429975603041186, 25.503197665951244, 26.15908554941181, 26.534137551883674, 27.6860326301588, 24.93114693726739, 27.848962647175245, 26.229235337926927, 26.40716029228892, 25.87154394254469, 23.65128117624798, 23.802465078586494, 24.23030948246349, 26.165343939774907, 25.034492788401177]),
array([27.77689612166788, 23.681778963985675, 23.21658161513045]),
array([26.21485216842258]),
array([25.408248013790836, 24.616032231669994]),
array([26.132483094305815, 26.17740924468248, 24.006147948903642, 22.875036703224612, 22.362786628505365, 22.45643381544724, 21.5383363918377, 21.583304416660866, 21.439657223211483, 21.12949102451826, 22.081611586520697, 21.314120986991245, 21.686127289338508, 22.907611630818074]),
array([24.32838708015632, 24.648007845187514, 21.9344141650573, 22.212343422290722, 23.005115246917747, 22.745147640328355, 22.52198093615798, 21.19341004213879, 22.861513588572084, 21.129532680129035, 21.700402122633037, 21.325448598493615, 20.562873068485697, 21.803531651423064, 22.08790124047093, 22.143683257001005, 22.331201334725925, 21.335656357647903, 17.0794099034719, 20.092756369828926, 17.92563223301481, 16.16867563992791, 17.24130767918431, 19.09045591792075, 17.92113419134647, 19.781826600479697, 16.504535051288435, 17.212061194665974, 17.083260789263843, 19.96332258703755, 16.12072366773017, 16.76497610888397, 19.67868143150915, 15.991701557076102, 17.66542562082534, 17.80970225504699, 20.314174034547154, 16.257901595180137, 20.43994127480161, 18.825622205678137, 17.057841480173217, 18.817380442440342, 19.81079824252609, 19.89703157169999, 16.560224150421885, 18.027732800120646, 16.124672561776652, 16.65717996092112, 19.572528553379467]),
array([24.196730052501366, 23.638204878956525, 23.05648697399193, 21.37134656918515, 22.60953591084721, 21.371204713169085, 21.06495271224985, 22.843286610466645, 21.814412845935852, 20.57706523934362, 21.681306314454236, 19.439401228185154, 19.071865121405477, 18.89728573456922, 20.265208755271654, 18.5432870716848, 20.075178504019135, 20.13671797775842, 19.56490033299593, 18.50188764723785, 18.58486792102086, 18.69274138149919, 18.761681579168897]),
array([23.61290615164405, 22.643159180765192]),
array([23.353120624276574, 22.65573338720521, 21.667319166342022, 22.175468913143114, 22.06567501281025, 16.587220633528883, 17.44212433302818, 19.28037059935815, 17.597791046684954, 20.376120117419735, 19.96171957522151, 15.786642167558384]),
array([20.887202196656347, 17.15960828980027, 20.034834469618282, 17.37613230717387, 19.623620800716807, 9.093722361442802]),
array([21.065943621272858, 21.60086039911353, 22.667521028473516, 22.153840994841186, 20.718539173416787, 20.014978767595228]),
array([21.296737859048957, 22.55899407781298, 22.240587048131253, 21.153750296404272]),
array([23.842324541442206, 23.74586859604844, 20.834917129837567, 21.29936315150979, 20.971268495845923, 22.82016106498413, 20.53976362070571, 22.011890119724242, 22.00274729709502, 22.2601033402082, 20.451329688729196, 20.769961423667162, 21.63314143531207, 21.084762631007035, 22.23831564683203, 18.851192046802296, 20.315355438742685, 16.75654684431779, 17.091187850522047, 18.796963928751964, 19.53318063700554, 16.22088253482108, 18.190996330281862, 16.634302973486314, 16.269932700081206, 16.170055176883444, 16.242860193253584, 19.77167336939638, 16.360499844110166, 18.473691626323976, 17.911463712008633, 18.49423579367251, 16.635558268659004, 19.147332871781547, 19.30615931935666, 17.571351335417276, 20.30982257359653, 19.89903408678608, 16.10588178611697, 19.704610294899442, 19.85279855468474, 18.8833468940341, 18.827831625340885, 16.01545398267721, 17.04965889566475, 16.18286628657873, 16.02490870414178, 16.453844400305478, 18.45405242646488, 16.971838226207833, 19.90935904316737, 18.400861903328863, 14.968574975734356, 15.09734379655599, 15.336282451977883, 15.713030067475476, 15.6066170080958, 15.935752225498323, 15.062062552125537, 14.110807016465166, 14.628969759697844, 12.630411662767665, 12.500048502282757, 12.382803124427, 12.00666522112927, 12.472393728618815, 11.176729665255996, 10.56987327351756]),
array([22.285943005574836, 22.51736350437002, 21.938064067502165, 22.49122802350679, 22.30267072649921, 21.995620406396903, 22.51965230882828, 22.946948324568535, 22.044466799482148, 22.30243898342781, 21.820532469251646, 22.517637563677, 22.661851473786744, 21.890786868334917, 22.777635304629126, 22.485407289608908]),
array([21.76225766575733, 20.992808515009177, 20.932986542368184, 21.685338318163375, 21.742301191918894, 20.507914810443012, 22.69715678709578, 21.650684714314473, 21.384221621120403, 21.039457076309816, 20.94447424245848, 20.353903140257998, 19.64782525720984, 18.954892367397573, 20.02141130428746, 19.486468445974094, 20.42323095352414, 19.30448562508589, 20.33179613709981, 19.599559830014705, 18.812461880325355, 18.915567147122132, 18.46041095444748, 20.215400458870583, 19.220512871976133, 18.427671240043477]),
array([21.017352731626154, 21.40645309271175, 21.8697860938547, 21.974463658897243, 22.465762450732587, 20.980266262057555, 21.0834248788019, 20.38810595622868, 18.119249003933046, 19.30756100733048, 17.58355163342663, 17.631597156850916, 19.615286189176373, 19.36021153557376, 20.009066774937875, 20.12869165121357, 16.66013242083092, 19.50737163696068, 20.298959744443053, 19.068894210539018, 17.511294812541358, 17.217022657939147, 16.649913135749177, 19.548176403026844, 20.438050976982638, 18.53954323061786, 19.224837481782554, 16.821933076735515, 17.5678339175333, 18.529591427309768, 20.405475395595406, 17.90904758409807, 17.64053049362113, 20.19676267507518, 18.69162675457202, 16.72302013850177, 17.207689882982756]),
array([20.948219596213548, 20.895481849758273, 21.493761378956524, 23.01093884639947, 22.905595286960047, 21.08058971727746, 19.19998647846832, 18.736178811747486, 19.951844048882215, 17.701471498903082, 18.16061070055769, 18.3681473770038, 17.748041171628834, 18.437222016043453, 19.047427133974832, 19.80173243231187, 20.323675204867648]),
array([23.430994101435537, 22.065910112300525, 22.815418152061405, 21.494685378591456, 22.225896282515425, 17.39233537660997, 17.258532993966224, 17.279357100825976, 19.266600121274017, 19.72770271629012]),
array([21.132603037823014, 22.10716319242673, 22.392506976143157, 21.286219799338408, 21.571305805316776, 22.241154161482594]),
array([22.89670745240747, 22.450601184522974, 21.625148190811423, 23.00228542894508]),
array([22.610413621192194, 22.963822896657767]),
array([22.759297093841443]),
array([20.998686076838933, 19.757039105157073, 19.319249438869573, 20.209858275248877, 19.781743160198605, 18.87157429252033]),
array([21.82309100452122, 21.828961043422684, 21.090884904170434, 21.22910936803755, 20.564011343170538, 21.344604458775557, 22.339021177296402, 20.84167285173541, 20.587765268177975, 20.69223366111454, 21.933057490384904, 20.6016849406014, 20.124364147891757, 20.252049036215855, 19.540428896244173, 19.94909673455338]),
array([22.139230914567083]),
array([20.46871590262925, 21.26092411615911, 20.795229085100964, 21.302132418499347, 21.496132308547455, 20.92638540531595, 21.682678944924213, 21.292531446816565, 21.566289087523174, 21.11422629285849, 21.28421263485823, 20.728895430315507, 21.150174505019475, 21.096007843056896, 20.517149828135143, 20.641959829972457, 20.625722213859223, 20.773848286262886, 20.55103623529288, 20.86726957595005, 21.45037645689199, 21.4728879548294, 17.069298195710918, 18.144315350770388, 18.547346169927675, 17.07319273769223, 16.69809372364822, 18.957690643293077, 16.50348722083744, 18.318712875211734, 16.616353055540593, 19.950199761241482, 18.084748015131744, 17.96580365575662, 20.124801502237993, 18.966026547807097, 17.5746332865818, 18.35763187567295, 16.862665104264217, 19.802445578317705, 20.245415965791008, 15.986748796813792, 16.316801599116015, 17.52581429348177, 19.34613809013339, 19.37681038573443, 18.576641269683176, 16.411177090912798, 19.474396695465764, 20.217321129478794, 19.671671161952858, 19.187557839564313, 16.708229975630328, 19.277193712108826, 18.901885876164062, 18.812133519779657, 18.989873486852623, 16.29146995632603, 17.919636375392486, 19.548625903955106, 19.171142683541394, 17.354296546896187, 20.072995435451162, 19.629227001632767, 18.85742486512852, 16.147718448352116, 16.573094663646977, 19.71584702827401, 16.894141582715513, 18.56070340306686, 17.46976807011985, 18.61744122419314, 19.15753086430653, 17.57576750412547, 17.358501209490907, 16.150973991533426, 19.146509200033414, 19.757135183483413, 17.425500308306333, 16.664503342241343, 20.32740204417148, 16.359694912139773, 19.6939322750276, 16.74827514979007, 17.774719610248614, 20.141620502611612, 17.66454818288009, 20.0739659974323, 18.62339277925044, 19.76594785253203, 20.416222017405772, 18.058288255361987, 20.432543975440158, 16.767257209986195, 20.009940629027884, 18.126040622004425, 16.4070566460857, 20.10933713099818, 16.0855759406802, 20.14141341074223, 19.277983744142738, 19.211013700742924, 18.52508034669099, 15.545299261030403, 15.545850100927863, 15.684892296042918]),
array([20.77690995821987, 21.015483327250443, 20.847663524220792, 20.776635662073723, 21.12787643615219, 16.35879026023241, 17.343368711017707, 19.80842424069713, 16.693580810548987, 18.31261966230788, 19.76755724122826, 19.86884438964567, 19.630625265753107, 16.682051732258223, 19.02029261558746, 17.492234232691775, 19.564400689461735, 20.173627767950364, 19.897064248738243, 19.1706543093972, 18.20741183574131, 17.507413581688812, 16.77387397038831, 18.745586248966298, 20.001591767048808, 18.889310483752705, 16.383846720838513, 15.997937360758126, 19.162935815035745, 16.091022435757214, 16.45912373376898, 18.858340058766647, 17.158597105582132, 19.959748530768728]),
array([16.92315614702504, 15.714016429216692]),
array([20.62534225830523, 19.75236607561958, 19.271912866931235, 18.642937697514682, 17.298209299227405, 17.325758205688444, 18.541504744654514, 19.705380201483162, 18.733269720832446, 17.898007100879372]),
array([19.383190812838713, 20.427988829824034, 19.347043231519226, 20.365047730031044, 18.17464415981186, 19.805147660459777, 18.62711409953998, 19.121230481876538, 18.560612289230754]),
array([16.48946315624413, 16.2602384474056, 16.662262104936243, 16.997794380027674, 17.007568447411643, 19.96622902066037, 17.83714266899771, 19.398924801393, 19.208295522542514, 17.033053281593684, 16.177805249617144, 18.352971190143137, 15.300227892221747]),
array([19.875434880511644, 19.404779328655167, 19.363819416222977]),
array([19.926961569417195]),
array([19.31415259131522, 18.69127711471894, 18.277240985873984]),
array([19.57919687536538]),
array([19.239228949197468, 17.842241098752595, 17.677693589750763, 18.236380166624524, 17.864112484184723, 18.788961262462536, 18.273519042289397, 17.997488114384556]),
array([17.212466493396594, 18.735262860443783, 5.323599456897857]),
array([16.888474930058013, 18.639909623478776, 18.37010470454806, 18.058480475978016, 19.076699772726172, 18.94159966104648, 18.18884598928347]),
array([18.86312696109468]),
array([18.4477209946333, 18.309284721615835, 18.420418175835366, 18.144608649492582, 18.894311327470422, 18.271201126389055, 18.52009480933107, 18.162883985769824]),
array([18.738515288982477, 18.57300742080703, 18.45566852698668, 18.96751201409776, 18.530217469404498]),
array([18.78026132749887, 18.06379768154113, 18.576125219513695, 18.776428816760053, 18.247401807546876, 18.63268562757985]),
array([18.459255286602986, 18.195332165251102, 17.979931327045445, 18.725075133968495, 17.803589100289386, 17.702647117720783, 17.76639993890695, 18.2248666317198]),
array([18.056298297331903, 17.483928861685282]),
array([18.07980118845218]),
array([17.865131703660307, 18.082854836908346, 18.31646152690633, 17.793508205093524, 17.961443650128587, 18.38894721721069]),
array([17.806786974826828, 16.849463441487053]),
array([17.349635060058645, 17.98636695443564, 18.01058552117745, 16.51673142315135, 16.93046349125338]),
array([17.19339303107028, 17.27291123428789, 16.968866342354232, 17.538925039135854, 17.171086786807038, 17.722522960338242]),
array([16.758482643061313, 17.42988934323196, 16.41342811328647, 16.34514898290075, 16.346059119429114, 17.274472285741144, 16.14218062427117, 16.22600383218249, 17.86154095891257, 17.863244867178572, 17.37767577351206, 17.69726594516461, 14.761123061061891, 14.837649064164056, 14.150475935402675, 15.57096959510196, 15.222819919017715, 15.2020312624316, 12.544376124377427, 13.05448910459224, 13.743324510844438, 13.781478115082184, 12.1018307537628, 12.34717867259496]),
array([17.631770397500958, 16.776625688554184, 17.21971827660359]),
array([17.353919308534465, 17.43173145601831]),
array([16.511221605636425, 16.577814368232573, 6.8733683111444]),
array([16.242096629597228, 16.26705414031547, 17.328471457316855, 16.343076603242462, 16.417195846581887, 14.505610295500539, 15.503681496039649, 15.965835993308058, 15.395330058069263, 14.516859510074692, 15.7586063772746, 13.80213041785537, 9.572599110621795, 9.276815760536524]),
array([16.497695313472207, 16.17231636139714, 16.575950167911802, 16.12311305250753, 17.333185908791563, 16.118085690316153, 15.789065222797435, 15.090578652687512, 14.827626205140847]),
array([16.194726070545276, 16.524549158487662, 16.97782935953082, 17.22425679033138, 16.855528627637113, 16.035999165197964, 16.922884442418656, 16.232599306907908, 17.26828963391641, 16.784359841764267, 17.26778798869021, 16.35146604006813, 16.321383685407273, 17.074842067220107, 17.325253474368314, 15.55512945090339, 15.28733662881195, 13.872277908997104, 15.943643218089244, 15.428380838597619, 15.681391193159678, 14.465827773673075, 15.811063195251105, 12.949731703376946, 12.769740180955958]),
array([17.039713937560162, 16.266298213117107, 16.44186996613531, 14.789642308327826, 15.222988822266176, 15.33387509387427, 15.498848936148764, 15.727987068781626, 15.611814709233784, 14.03786668167256, 11.800284186675093, 12.859452836750815, 12.059220840829783, 12.244707003325544]),
array([14.683375515958312, 15.74728037321395, 8.483646118908432]),
array([16.902074262768707]),
array([16.0154093159945, 16.616484033590513, 16.130415849255286, 16.695240091256277, 16.010428600787964, 16.207794908155336, 16.437297624239257, 16.626675988573844, 16.792019797800904, 16.351492233213367, 16.361534302754105, 15.406646088463429, 15.151644514357464, 15.50375608197488, 15.19667871205709, 15.858080190250089, 15.614958785762013, 15.385130966138455, 15.753290513732466, 15.080272452213459]),
array([16.623569504869014, 16.774379323963217, 16.846270458661, 16.807230363837178]),
array([16.57488593843995, 16.241471396624306, 14.178443078802113, 14.657261984678454, 13.896280956861823, 15.029506133710068, 14.445653461484486, 15.414804689429566, 14.598875874252208, 13.60138594575035]),
array([16.783586756336224, 16.801531729453608, 16.856628853776485, 16.787624293013845]),
array([16.149379075246053]),
array([16.050517747004854, 16.762080141882045, 16.342289581838724, 16.326858604461712, 16.084503914724355, 16.77032118484373, 16.002452252485117, 15.222906035516393, 15.897403273768251, 15.30873932328325]),
array([14.870353692180277, 14.84856298680167, 14.209060057787118, 13.973025277062277, 15.726491464348115]),
array([14.643647589105033, 15.75222015536091, 15.384753344772015, 14.435201281416536, 15.057596818431035, 15.213024484320469, 15.034333343180629, 14.473357199225678, 13.120088873720606]),
array([15.97784452243432]),
array([13.655922882991506]),
array([14.68965264916783]),
array([14.50686105683275, 13.864085973285212, 13.823027144903788, 14.386626373307667]),
array([13.838453085295235, 13.40809879365727, 12.541161317671904]),
array([11.616179426218256, 7.3475290443117265, 11.624872876735214, 5.165785611038676]),
array([13.171670708068776, 11.113365637299335]),
array([12.650830130009982, 12.865810645828532]),
array([11.656883736015738, 12.577452721951396, 12.323409899238301, 12.125268565496159, 11.955464349567782, 11.18056092040374, 9.668135738673824, 11.387570926647774, 11.337008756962334]),
array([12.47811701145651, 12.790954451440083, 12.521082357495928]),
array([11.846476485506205, 12.648153339056346, 11.791885528173394, 11.82536524141228, 11.788907809002069, 10.450228093037603, 11.301867943375928, 10.910656908086699, 10.429259072952485, 11.504055446331837, 10.040291738355018, 11.1019050767713]),
array([10.516256412857416, 9.17491074356142, 7.128855341429062, 4.761624730188029, 2.7263567943275944, 2.7563872364216575, 1.5826737577329963]),
array([11.085825885245063, 11.4129149698238, 7.64249016707244, 8.139536296143127, 7.3255935175598905, 8.553755607069656, 7.171094869361423, 6.890184860412378, 5.973360002257296, 5.809127437817205, 5.533175451370763, 6.621203968968735, 5.869409669415255, 5.661218006055595, 6.172968796042338, 7.123451505846034, 3.7741701166864994, 5.164732335057343, 4.877984911755395, 3.735675616479853, 4.488022345604375, 4.981134669143643, 3.800357464763472, 3.2340896329603193, 3.456471928446617, 3.1292995139678466, 2.7163395786072435, 2.586695660845371, 2.666594192731917, 3.364627362326148, 2.827260714033568, 3.422033330722079, 3.109068670438345, 3.3431269074277985, 2.92235211728464, 3.251022075653239, 3.2511069475984984, 2.634055677347191, 2.4983459124043708, 2.169615590341552, 2.456187957848728, 2.051745785359431, 2.479097782637931, 1.5153392751041406, 1.3167663631011186, 1.494357046690055]),
array([7.824942465834608, 11.541342993757105, 6.916719936397648, 3.055146481823835, 2.972760240311831, 3.023759807998788, 2.9336791853580397, 2.931046132089496, 3.303800937636666, 3.3671346954150776, 3.5967543584427153, 2.4228100037513194, 1.562681886787381, 1.2071687975526764, 0.3145450602900591, 0.7148968448247287, 0.4412763215745703, 0.0]),
array([11.757334267378381, 11.337910973635791, 10.502395434916796, 10.73819985141594]),
array([11.76635466794904]),
array([9.357047149823012, 8.322710379914017, 9.065188077326582]),
array([9.724978561851763, 0.0]),
array([10.85817110364878]),
array([10.356768165503956]),
array([10.66870676598393]),
array([9.882854879052854, 8.89187726112992, 10.489023892791053, 8.71156776264929, 7.976995558862433, 10.160427809783291, 8.179889567321446, 10.460370371340582, 10.787024071140785, 5.985137724435965, 6.500186502148613, 6.351212954629337, 5.490201752284286, 5.610041716766995, 6.739898190397618, 5.491913929282856, 6.131228245174892, 6.370496340971381, 7.086656211728565, 6.054208377048056, 5.530555529635307, 6.238996034711901, 6.178060461949682, 5.1635453092059125, 5.151325337800813, 4.955748728706956, 4.22453189469591, 4.413845832235949, 4.206773018101555, 5.225312250663396, 4.705280855335005, 2.587741085136577, 2.950391277724835, 3.0374026315942535, 2.7514707416577955, 3.3097720117568694, 3.2777609741387668, 3.0300671187666115, 3.5875248652854554, 2.6530722932460833, 2.784853023405025, 3.536557757113641, 2.8315027988419645, 2.9945694780597254, 2.415390511920762, 2.3337777882167807, 2.492891533782061]),
array([7.67128333280124, 9.690672118354842, 7.3491335980425605, 9.250272127823804, 7.021561912516951, 6.864864158705981, 5.477242976715745, 6.10702168822764, 6.950241119048458, 5.421049449476652, 6.158569560137701, 3.9831755348258264, 4.307191858205116, 2.769136669271923, 2.7119265965049113, 3.2320986653988326, 3.572118850297407, 2.9420465906992437, 3.1342200025768516, 2.3719689017592427]),
array([10.396670002421057, 9.223475638624233]),
array([9.911322393860376]),
array([8.900594753886308, 9.155592027597237, 8.93876781304817, 9.824829338244173]),
array([10.232181372392384, 8.190490264089231, 6.879118172178123, 7.151283906875355]),
array([8.099563131190598, 9.167751652804691, 7.5641054074553935, 7.778815765402062, 6.556770727880513, 6.130562821255484, 7.11323397746898, 6.723785181641666, 3.656471651106238, 4.356213301443091, 4.729317915285997, 3.9777697217525736, 4.486915417194018, 3.528160921845769, 3.512223743699991, 2.8690807196583474, 2.6968359477618797, 3.1900337131011707, 2.828731066528171, 3.4513347198800512, 3.4615557815333875, 3.458956883359619, 2.013280092515762, 2.434175744766835, 1.6389970116617938, 1.5052172023645436, 1.6206434557662641, 1.139249725205718, 1.0637905046849434, 0.9393129076325463, 1.6998426761706136, 0.631524775073341, 0.2740402134696739, 0.0]),
array([7.347208881569997, 9.510912551163017, 6.858424052161171, 6.6534310921793365]),
array([9.218921672112568, 5.386290007931356, 7.227077750992431, 6.044640058464816, 6.843348918604538, 5.426405979943914, 6.632591003205147, 5.715176825469042, 7.09409176969069, 2.906423018376048, 2.8684632828117773, 3.1771216124877606]),
array([9.114743389501397, 6.87827334710403, 6.09162394809811, 2.7491438187018344, 3.228128684180102, 2.912942738361819, 3.044738347910134, 2.4402418829576225, 1.9211499946066621, 1.8735884678088417, 1.6404156614323337, 1.42591400101959, 1.2643083593975948, 1.1772017141375604, 0.6732831916018991, 0.6806442980213068, 0.19437471070328827, 0.0]),
array([8.230213500719756, 5.669152972167971, 6.038014452818214, 4.792517299580818]),
array([9.008583175739268]),
array([7.826124696721086, 7.377179150620507, 5.885105501210592, 6.289183398586589, 7.236623612207634, 7.165280366095421, 6.773349541477681]),
array([7.500588356818484, 5.3857543128259575, 7.2101704017862245, 3.900259106178653, 4.7908678737669295, 4.376735397524553, 3.4651037470471633, 3.289456877412055, 2.8602016293300974, 2.523997951970402, 2.444553202260296, 1.4718857175952331, 0.5555236757912333, 0.6858425244320864, 0.31851629091151407, 0.26782287864530263, 0.4791095494932966, 0.5406979724796792, 0.0]),
array([7.508484400584821, 8.24543087390224]),
array([7.277944310535102, 6.65754949116222, 5.34964567470241, 6.059302887365163, 5.853778696890151, 5.757398136780246, 5.128374235104559, 4.9109047230200495, 3.569698088696712, 2.696950859751051, 1.517071143073685, 1.372273836980067, 0.3819597346990299, 0.37921808168405374, 0.0018699058084883785, 0.0]),
array([3.8164166393261736]),
array([7.864877073435876]),
array([7.673333030208724, 7.394789215317642, 6.862472671372694, 5.931393773149404]),
array([5.631310661989296]),
array([6.884088082629648, 2.593483317664015, 2.9572228168001735, 1.8213612752062747, 2.341043799419257, 1.530802536221017, 0.0]),
array([7.404920156200827, 5.82328246187325, 6.867385155246353, 6.958433789713998, 7.166282999955728, 5.454268798862513, 6.5548812105254255, 7.245538691965719, 6.936093269577205, 4.862408373710326]),
array([1.4494910997957378, 0.0]),
array([5.362148474539492, 2.23724376735035]),
array([7.350671787981774, 7.293405861351666, 7.143476410319232, 7.029163631514846]),
array([7.053002014779125]),
array([7.100013772938965]),
array([7.159874156163508, 6.143198330746535, 7.009736882777218, 6.715439049308368, 5.335683014109032, 7.141826597098223, 6.39785799345057, 5.636983709068943, 6.844526896901296, 5.555868582538347, 6.86967703487139, 6.008661151227232, 5.311127571959276]),
array([5.768022770319389, 4.119130391587931, 4.309926316358057, 3.3568305504819174, 1.9548899280386167, 0.3624526271326206, 0.6648426052744892, 0.0]),
array([6.3898858755178605, 5.70959338393734, 6.5982145412693605]),
array([3.5138843386069487]),
array([6.489477829586328]),
array([5.134050304783789, 2.218580892036134, 1.927232504450103, 2.031370041257697, 2.2031600159490186, 1.1298874656306392, 1.2834858189750373, 0.5490325099767729]),
array([5.382290505541513, 2.65351745994598, 3.13900604880264, 2.4003474445455373, 0.8649529507523535, 0.23028105076444694, 0.03371417327044671, 0.0]),
array([5.8661732199975996, 3.319803647458816, 2.6388953108763986, 2.670024890517994, 3.146038297244732, 2.8548694453791326, 1.3407200409572613, 0.7680643027425945, 0.6818192473722158, 0.0]),
array([6.019347723247378, 3.547659440571773, 2.6876745741497396, 3.2532246125300746, 3.0126058265969333, 1.9029699720989053, 2.289560326425586, 1.2770043959505277, 1.6931665879354307]),
array([5.730907889090596, 5.77940353772442, 4.862407036272942, 3.860963349314093, 3.578607216303614, 3.04213168246958, 2.367985959125715, 2.5007040679312, 1.8833057614919144, 1.637369982666681, 0.8806571168389096, 1.708972154498216, 1.665820572564606, 0.6274993139896645, 0.14020181302061807, 0.0]),
array([3.5801107688645786]),
array([5.626549341748123, 2.6996804924718862]),
array([5.388431700544497]),
array([5.375271501316051, 5.422494803971996, 5.154752680467319, 5.3272999047244, 4.457842151530581]),
array([3.31813786507417]),
array([5.331487075296265, 4.965011775985767]),
array([3.730804291220119, 4.861902229803265, 3.4555360968513558, 3.485289518974357, 3.520100849267662, 2.8505333895793443, 2.9072296481258517, 3.1699527482124608, 3.5987001586475578, 3.137266696322043, 1.8674938336020768, 2.4731231430045852, 2.1162975026575905, 1.5468394611632479, 1.6818393335161925, 1.5488034988152923, 0.8724715915715029, 1.6144212675513494, 1.7073373798456388, 1.6288882428807379, 0.321794468951141, 0.1605792782543296, 0.5514071005175522, 0.48506112275794944, 0.31474062460550556, 0.4229449054822318, 0.3744741351103797, 0.15246678688648962, 0.05643270222484795, 0.0]),
array([4.452929433505668, 4.317685280951938, 2.6120803010817006, 2.8540334850995825, 1.6278539612759462]),
array([4.93588609183684]),
array([3.144211721136716, 3.02433911556926, 0.8456226759038981, 1.3836591185322804, 1.6337934409210428, 1.0761927051655142, 1.3975932256002013, 1.186764705136544, 0.0]),
array([2.803895681621185]),
array([2.802051472081347, 1.847530812800323]),
array([3.290593043930305, 2.706118991152739, 2.604932369360189, 2.760424511935104, 2.8153413190640086, 2.1207965133119155, 1.9890687325258867, 1.296700617738162, 1.3864129692873224]),
array([4.387802322966189, 2.8310043823383637, 1.8187912000136666, 2.332128821094059, 0.0]),
array([4.043400940135573, 3.54538212853258, 1.7334819375507053]),
array([1.691227962031168, 0.2766265035467057, 0.3441133023397087, 0.12636990749523214, 0.0]),
array([4.020519387526392, 3.470661982654633]),
array([4.401912563997259, 4.469076743797847, 4.510636638817477, 4.510340772089646, 4.199377734786373, 3.9568989702511255, 3.9402642755488735, 2.6666707726814463, 3.027537287565866, 2.9097742766685895, 2.7270508968442035, 3.178544826367905, 2.5859406073028604, 3.077850305700212, 1.8994129948749556, 2.1603406258474562, 2.5382393105673597, 1.5305137811079685, 1.6082509881756533, 1.641649330556414, 1.364887973077926, 1.6723696096733847]),
array([4.02546644601054, 3.8535264926100794, 2.7877697474062395, 3.361918357336163, 2.983123016581663, 2.870042174520038, 2.738446939257752, 2.4059282632579024, 0.15377007332780346, 0.0]),
array([4.385935055911334, 3.799819576879174, 4.043272571143368, 3.488348253972193, 3.2843049225058825, 2.9594906775599465, 3.1173327453429325]),
array([3.682826713613891, 3.1882465591772395, 3.5225982819720203, 3.5136858057027913, 3.5277223644031728]),
array([3.1692094021619615, 2.91856460984455, 2.4029767134770093]),
array([4.227029640395331, 3.89808920264614, 2.790951773723986, 3.0854489655204334, 2.6803713110655893, 3.3521515910735995, 2.959338977620101]),
array([3.4584209436568085, 3.5203260114107486, 3.3973682131413816, 3.1528114212662546, 3.0350051341782325, 2.6708352867669296, 3.3027460281666294, 2.3788839513480067, 2.461064408875212, 1.1555895896966644, 0.8057668776979087, 1.331187634552068, 0.9694230958038075, 1.468166935065215]),
array([3.2625632444590793, 2.6315412180044535, 3.530851316104599, 2.309172154514897, 1.9233818082337062, 1.9387330035389136, 1.460105411989005, 0.586006944620497, 0.0]),
array([3.5955254303281654]),
array([1.8844693622290714, 2.260416379811372, 1.4185661979857407, 0.19784722553185063, 0.0]),
array([3.2395172042987794, 2.661342302999084, 2.8247641545463007, 2.812693332619509, 2.4992336068792915, 2.3381433633731734, 1.8522481994073647, 2.542853079880854, 0.9795072975890438, 1.4987705173682593, 1.0886133559695041, 0.7152109375364686, 0.6188804089454416]),
array([2.8837239358180677, 3.319295164769035]),
array([3.083520567902143, 2.6421050463641267, 3.1429548453450007, 3.120023804576859, 2.4828220279314195, 2.1797817394186874, 2.404068891024555, 2.1189054298956496, 1.96123242146877, 1.507434600149326, 1.51587208438526, 1.6863740658480404]),
array([2.65577922899275, 2.124298414840861, 1.8673407460285905, 2.1625307277429497]),
array([2.643518953794467, 2.9432296863731455, 2.970889944763794, 2.868131306530969, 3.006978179356637, 2.622518272877281]),
array([2.936864818060257, 2.618577576738411, 2.928948581105541, 2.57009344113391, 0.7763661248408932, 0.6081594666937393, 0.0]),
array([2.588456934771623]),
array([2.658543318444238, 2.667898268725253, 2.502796331830433, 1.5681992352940366, 0.9142072514925464, 1.4345455336930075, 1.6668348601664988, 1.127052613666612, 1.0615824571080117, 1.0837861311530492, 0.2281497750892738, 0.3228334242802158, 0.0]),
array([2.543027471209289]),
array([1.936225753598415, 2.0611518912061415, 2.3324389784699013, 1.8795605471329264, 2.379558026366931]),
array([0.13451655814423502, 0.45180281619580287, 0.0]),
array([2.062850270170144, 0.4107534997783897, 0.43318182312335846, 0.0]),
array([2.0882466760611167, 1.070066578534437, 1.6165022953490702, 1.4953076824112195, 1.669478391507935, 1.6619387529925438, 0.6436783808723912, 0.45100598775721396, 0.0]),
array([1.8320249244894578, 1.6962667108139942]),
array([1.1566351313909378, 1.0416451006530898, 1.2429034072199263, 1.0255397733441236, 1.676574022108263, 0.6126626065392584, 0.5716570302489358, 0.2866725801251316, 0.68325533302319, 0.4505571777813737, 0.0]),
array([1.721824377609285, 1.4430997022580052, 1.2732587188118583, 1.617920914924828, 0.7687830273109141, 0.6596435496440277]),
array([1.0562987809192759, 0.9695851604495797, 1.4875482834485896, 0.4812771316922391, 0.16039651127368426, 0.446420312099622, 0.005417742712934881, 0.048053908954367466, 0.12456140499649736, 0.0]),
array([1.0777661288965432, 1.2897736594262275, 0.47269936313008415, 0.0]),
array([0.9316022056986888, 1.1832421806101188, 1.0459753991589267, 1.3488556242638796, 1.348845990157626, 0.9133786947017746, 1.1292507803569838, 0.794098574951982, 0.16390551985442448, 0.31377332514658407, 0.23042031296245247, 0.22918546358804714, 0.6352913097732427, 0.053509598710937134, 0.0]),
array([1.552599641888752, 1.5277042725319385, 1.1828126511779888, 0.7873251192596652]),
array([1.1563607118702868, 1.1620444128598824, 1.1935518702838395, 1.2726703794013945, 1.3797700129095283, 1.118928323037503, 1.1342740383644034, 0.5438727060533446, 0.33433622593079143, 0.38357423094222887, 0.7651968313537939, 0.24281984089700348, 0.32964315850276327, 0.5085803573079359, 0.5241713157167691, 0.3514441165584664, 0.7123268600583743, 0.7083898181554199, 0.4515093439236829, 0.3294774733491307, 0.4924133814487043, 0.6446330688216073, 0.7516474191469613, 0.15967513547452816, 0.10002614807259445, 0.0]),
array([1.0921953309368622, 1.1645562779242895, 1.0040608707379521, 1.1553278451825395, 1.2558566804785605, 1.1287541105004661, 1.0577976302802274, 0.1770112897387205, 0.578154876153992, 0.3503361363523106, 0.5584144531791868, 0.485641310558414, 0.5777473334795094, 0.7131164527402643, 0.4104606765956916, 0.1591984006896301, 0.217491724094095, 0.0]),
array([1.097352389523664, 0.9809040644487517, 1.015854239258353, 1.1457702805364602, 0.7942967464142774, 1.1657587157488778, 1.0483357079286906, 1.0450818338483845, 1.2527542470134725, 0.7026705239652634, 0.2615321949466807, 0.33428069352622375, 0.12883677063488697, 0.171447627340652, 0.5930905898520922, 0.36146114715768385, 0.19228242998547107, 0.30395650872693214, 0.4326247534087936, 0.6367545846383965, 0.39406200963795873, 0.3606806151535429, 0.4761542931760827, 0.2643250055908235, 0.7402869288760424, 0.5501606279353927, 0.7796072208420303, 0.21066211501249632, 0.3998086460148324, 0.68993666390084, 0.48825656304340953, 0.4399940334310535, 0.49045226896606836, 0.0]),
array([1.0943347972712651, 0.3393798066970124, 0.41382953015184926, 0.4258115587334459, 0.35405345918762143, 0.0]),
array([0.31189899480321703, 0.3208516236335359, 0.7362683465787243, 0.0]),
array([0.22558463276011453, 0.7250853882256298, 0.0]),
array([0.9728249893203782, 0.6235318292878577, 0.3398572985175868, 0.43202237046108255, 0.31705515720995825, 0.19283351040323637, 0.0]),
array([0.8724479754266431, 0.4634597053998869, 0.3493902559530879, 0.5175665405405584, 0.4689575509802288, 0.3555518138645531, 0.4565041405076198, 0.7014053444324696, 0.6590075295034937, 0.0]),
array([0.6435162056022731]),
array([0.8168924933875182, 0.7101514203853473]),
array([0.6855908592715653, 0.2598955271478487, 0.5682228029630845, 0.6137184168679815, 0.0]),
array([0.5736881031328183]),
array([0.5455276356709078, 0.21725825052095704, 0.6003075481703305, 0.4634442442890134, 0.2543131193941107, 0.0]),
array([0.329152380205001, 0.253431411784345, 0.48674529028253755, 0.5769993874104112, 0.4663886662379733, 0.0]),
array([0.09307664753989867, 0.0]),
array([0.2360832262546624, 0.19286664856677316, 0.26099016891762505, 0.0]),
array([0.13955624503068145, 0.1168104016998678, 0.0])
]
d = [data_1]
names = ["75"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T2', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T11', 'T12', 'T13', 'T15', 'T16', 'T17', 'T19', 'T20', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T34', 'T37', 'T38', 'T40', 'T41', 'T43', 'T45', 'T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T86', 'T87', 'T88', 'T89', 'T90', 'T92', 'T93', 'T94', 'T95', 'T99', 'T100', 'T101', 'T103', 'T106', 'T107', 'T108', 'T109', 'T110', 'T112', 'T113', 'T114', 'T115', 'T116', 'T118', 'T119', 'T121', 'T122', 'T125', 'T126', 'T127', 'T131', 'T132', 'T133', 'T135', 'T136', 'T138', 'T139', 'T140', 'T141', 'T145', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T157', 'T158', 'T161', 'T163', 'T164', 'T167', 'T171', 'T172', 'T173', 'T175', 'T176', 'T177', 'T178', 'T179', 'T181', 'T182', 'T184', 'T186', 'T187', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T198', 'T199', 'T200', 'T202', 'T203', 'T204', 'T207', 'T208', 'T209', 'T212', 'T213', 'T214', 'T215', 'T216', 'T217', 'T219', 'T221', 'T222', 'T223', 'T226', 'T227', 'T228', 'T229', 'T230', 'T231', 'T232', 'T235', 'T237', 'T238', 'T239', 'T240', 'T241', 'T242', 'T245', 'T246', 'T247', 'T248', 'T249', 'T250', 'T252', 'T254', 'T255', 'T257', 'T260', 'T262', 'T264', 'T265', 'T266', 'T268', 'T271']
def get_taxa_names(): return taxa_names