function retVals=BernoulliNum(n)
%%BERNOULLINUM Obtain the nth Bernoulli number. Bernoulli numbers play a
%           role in Bernoulli polynomials, which play a role in the
%           development of lattice rules for numerical integration.
%           Bernoulli numbers also appear in formula for the explixit
%           expression of the sum of powers of the first n positive
%           integers. The definition for BernoulliNum(1) used here is -1/2;
%           some authors use the opposite sign.
%
%INPUT: n The index of the desired Bernoulli number. This can be an array
%         or matrix if multiple indices are desired. The indexation starts
%         from 0, as is the traditional subscript used.
%
%OUTPUTS: retVals The Bernoulli number(s) at the positions given in n.
%
%Numerous algorithms for computing Bernoulli numbers exist. However, most
%of them are problematic, due to finite precision limitations. Thus, the
%best way to compute such numbers is to use a symbolic algeabra program and
%simplify the results. However, IEEE double-precision floating point
%numbers can only represent even values of Bernoulli numbers up to 258 (odd
%values > 1 are 0). Thus, this function just reads the values from a table
%that was computed using an arbitrary precision library. Every value in the
%table is represented as an exact fraction.
%
%An expression for the 2n-th Bernoulli number that is good for use with
%algorithms for symbolic solvers is Algorithm 2 in [1]. As noted ion the
%paper -2*nzeta(2*n+1) given the 2*n-th Bernoulli number. Algorithm 2 in
%the paper when run with 2*n terms provides the exact solution.

%REFERENCE:
%[1] P. Borwein, "An efficient algorithm for the Riemann zeta function," in
%    Canadian Mathematical Society Conference Proceedings, vol. 27, 2000.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The table of all Bernoullli numbers in rational form up to where the size
%exceeds the magnitude of the largest IEEE double floating point number.
valTable=[1;%n=0
1/6;%n=2
-(1/30);%n=4
1/42;%n=6...
-(1/30);
5/66;
-(691/2730);
7/6;
-(3617/510);
43867/798;
-(174611/330);
854513/138;
-(236364091/2730);
8553103/6;
-(23749461029/870);
8615841276005/14322;
-(7709321041217/510);
2577687858367/6;
-(26315271553053477373/1919190);
2929993913841559/6;
-(261082718496449122051/13530);
1520097643918070802691/1806;
-(27833269579301024235023/690);
596451111593912163277961/282;
-(5609403368997817686249127547/46410);
495057205241079648212477525/66;
-(801165718135489957347924991853/1590);
29149963634884862421418123812691/798;
-(2479392929313226753685415739663229/870);
84483613348880041862046775994036021/354;
-(1215233140483755572040304994079820246041491/56786730);
12300585434086858541953039857403386151/6;
-(106783830147866529886385444979142647942017/510);
1472600022126335654051619428551932342241899101/64722;
-(78773130858718728141909149208474606244347001/30);
1505381347333367003803076567377857208511438160235/4686;
-(5827954961669944110438277244641067365282488301844260429/140100870);
34152417289221168014330073731472635186688307783087/6;
-(24655088825935372707687196040585199904365267828865801/30);
414846365575400828295179035549542073492199375372400483487/3318;
-(4603784299479457646935574969019046849794257872751288919656867/230010);
1677014149185145836823154509786269900207736027570253414881613/498;
-(2024576195935290360231131160111731009989917391198090877281083932477/3404310);
660714619417678653573847847426261496277830686653388931761996983/6;
-(1311426488674017507995511424019311843345750275572028644296919890574047/61410);
1179057279021082799884123351249215083775254949669647116231545215727922535/272118;
-(1295585948207537527989427828538576749659341483719435143023316326829946247/1410);
1220813806579744469607301679413201203958508415202696621436215105284649447/6;
-(211600449597266513097597728109824233673043954389060234150638733420050668349987259/4501770);
67908260672905495624051117546403605607342195728504487509073961249992947058239/6;
-(94598037819122125295227433069493721872702841533066936133385696204311395415197247711/33330);
3204019410860907078243020782116241775491817197152717450679002501086861530836678158791/4326;
-(319533631363830011287103352796174274671189606078272738327103470162849568365549721224053/1590);
36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901/642;
-(3469342247847828789552088659323852541399766785760491146870005891371501266319724897592306597338057/209191710);
7645992940484742892248134246724347500528752413412307906683593870759797606269585779977930217515/1518;
-(2650879602155099713352597214685162014443151499192509896451788427680966756514875515366781203552600109/1671270);
21737832319369163333310761086652991475721156679090831360806110114933605484234593650904188618562649/42;
-(309553916571842976912513458033841416869004128064329844245504045721008957524571968271388199595754752259/1770);
366963119969713111534947151585585006684606361080699204301059440676414485045806461889371776354517095799/6;
-(51507486535079109061843996857849983274095170353262675213092869167199297474922985358811329367077682677803282070131/2328255930);
49633666079262581912532637475990757438722790311060139770309311793150683214100431329033113678098037968564431/6;
-(95876775334247128750774903107542444620578830013297336819553512729358593354435944413631943610268472689094609001/30);
5556330281949274850616324408918951380525567307126747246796782304333594286400508981287241419934529638692081513802696639/4357878;
-(267754707742548082886954405585282394779291459592551740629978686063357792734863530145362663093519862048495908453718017/510);
1928215175136130915645299522271596435307611010164728458783733020528548622403504078595174411693893882739334735142562418015/8646;
-(410951945846993378209020486523571938123258077870477502433469747962650070754704863812646392801863686694106805747335370312946831/4206930);
264590171870717725633635737248879015151254525593168688411918554840667765591690540727987316391252434348664694639349484190167/6;
-(84290226343367405131287578060366193649336612397547435767189206912230442242628212786558235455817749737691517685781164837036649737/4110);
2694866548990880936043851683724113040849078494664282483862150893060478501559546243423633375693325757795709438325907154973590288136429/274386;
-(3289490986435898803930699548851884006880537476931130981307467085162504802973618096693859598125274741604181467826651144393874696601946049/679470);
14731853280888589565870080442453214239804217023990642676194878997407546061581643106569966189211748270209483494554402556608073385149191/6;
-(3050244698373607565035155836901726357405007104256566761884191852434851033744761276392695669329626855965183503295793517411526056244431024612640493/2381714790);
4120570026280114871526113315907864026165545608808541153973817680034790262683524284855810008621905238290240143481403022987037271683989824863/6;
-(1691737145614018979865561095112166189607682852147301400816480675916957871178648433284821493606361235973346584667336181793937950344828557898347149/4470);
463365579389162741443284425811806264982233725425295799852299807325379315501572305760030594769688296308375193913787703707693010224101613904227979066275/2162622;
-(3737018141155108502105892888491282165837489531488932951768507127182409731328472084456653639812530140212355374618917309552824925858430886313795805601/30);
10259718682038021051027794238379184461025738652460569233992776489750881337506863808448685054322627708245455888249006715516690124228801409697850408284121/138;
-(81718086083262628510756459753673452313595710396116467582152090596092548699138346942995509488284650803976836337164670494733866559829768848363506624334818961419869/1794590070);
171672676901153210072183083506103395137513922274029564150500135265308148197358551999205867870374013289728260984269623579880772408522396975250682773558018919/6;
-(4240860794203310376065563492361156949989398087086373214710625778458441940477839981850928830420029285687066701804645453159767402961229305942765784122421197736180867/230010);
1584451495144416428390934243279426140836596476080786316960222380784239380974799880364363647978168634590418215854419793716549388865905348534375629928732008786233507729/130074;
-(20538064609143216265571979586692646837805331023148645068133372383930344948316600591203926388540940814833173322793804325084945094828524860626092013547281335356200073083/2490);
5734032969370860921631095311392645731505222358555208498573088911303001784652122964703205752709194193095246308611264121678834250704468082648313788124754168671815815821441/1002;
-(13844828515176396081238346585063517228531109156984345249260453934317772754836791258987516540324983611569758649525983347408589045734176589270143058509026392246407576578281097477/3404310);
195334207626637530414976779238462234481410337350988427215139995707346979124686918267688171536352650572535330369818176979951931477427594872783018749894699157917782460035894085/66;
-(11443702211333328447187179942991846613008046506032421731755258148665287832264931024781365962633301701773088470841621804328201008020129996955549467573217659587609679405537739509973/5190);
4166161554662042831884959593250717297395614318182561412048180684077407803317591270831194619293832107482426945655143357909807251852859279483176373435697607639883085093246499347128331/2478;
-(1369347910486705707645621362512824332220360774476594348356938715366608044588614657557436131706543948464159947970464346070253278291989696390096800799614617317655510118710460076077638883999/1043970);
1124251816617941290026484851206299982774720467712867275292043701618829826708395745459654170718363182143418314514085426692857018428614935412736063946853033094328968069656979232446257101741/1074;
-(6173136454016248924640522272263470960199559328290655337530202055853397791747341312347030141906500993752700612233695954532816018207721731818225290076670213481102834647254685911917265818955932383093313/7225713885390);
4277269279349192541137304400628629348327468135828402291661683018622451659989595510712915810436238721139546963558655260384328988773219688091443529626531335687951612545946030357929306651006711/6;
-(857321333523056180131194437347933216431403305730705359015465649285681432317514010686029079324479659634642384809061711319481020030715989009140595170556956196762318625529645723516532076273012244047/1410);
22258646098436968050639602221816385181596567918515338169946670500599612225742487595012775838387331550474751212260636163500086787417640903770807353228157478339547041472679880890292167353534100797481/42;
-(14158277750623758793309386870401397333112823632717478051426522029712001260747920789473711562165031101665618225654329210473605281619696918061316240634857984019071572591940586875558943580878119388321001/30);
5411555842544259796131885546196787277987837486638756184149141588783989774511509608733429067517383750706299486822702171672522203106730993581242777825864203487238429479957280273093904025319950569633979493395/12606;
-(346465752997582699690191405750952366871923192340955593486485715370392154894102000406980162521728492501917598012711402163530166516991115122131398542029056286959857727373568402417020319761912636411646719477318166587/868841610);
2269186825161532962833665086968359967389321429297588337232986752409765414223476696863199759981611817660735753831323900456495253961837175924312108872915089534970310604331636484174526399721365966337809334021247/6;
-(62753135110461193672553106699893713603153054153311895305590639107017824640241378480484625554578576142115835788960865534532214560982925549798683762705231316611716668749347221458005671217067357943416524984438771831113/171390);
88527914861348004968400581010530565220544526400339548429439843908721196349579494069282285662653465989920237253162555666526385826449862863083834096823053048072002986184254693991336699593468906111158296442729034119206322233/244713882;
-(498384049428333414764928632140399662108495887457206674968055822617263669621523687568865802302210999132601412697613279391058654527145340515840099290478026350382802884371712359337984274122861159800280019110197888555893671151/1366530);
2250525326187264545900714460628885135841050444551247116222631411681549780530233516069957534394574922579290608180427520318235621123686109474343887857944611842438698399885295153935574958275021715116120056995036417537079471/6;
-(110636644250856903590976481422794879200517231299540994715372334521128669716264196333811025709747746193210786820114369025849897345722531098042760530922656878891556664782168465095563132092311332073097630676251482491663634626858373/281190);
2525292668891404920279427026668969389456388249389889339455604316691573384284678293620100066924361693666444722338743839198221347931651916807651198800935942493038194104759967208073711284671045255047521429204396148980705984836743/6;
-(12407390668433023412711473483696990726334795896412761472587854072142800403373577087021298541061094633377354326966623278849423631924808044397822651135905640812063181221280972334965193338438214107578486417026806166184210160001817890901/27030);
4708181368529492614110644197951837317202610608341257204206693195241245204360822875910613010433572133227831741097261618833206537519198885812254347219150482005543422997225440204041473518187636442241332621804718967775203938403965710395632762155/9225988926;
-(1856110669947388268389361040689764027464160460436671923253131176853224087741924378432403442710398247642246902212818749685974336641529240178398124235555437625251481044526024910356131819016670047949661636539964662370375622630863327168696307/3210);
4005748930070152861935826766476856180706477227448622268042052745245798242539770546339789899546160341590069109467023517085578618986055969187202731878271685432460708841118506310943865592568791360294244451765746911808994482063783730693703607/6;
-(11993122770108617858536443322964878003618156069559794803117809279608039120818829088000103355036592864877954563564831932363414886811786054601318517206937549605059298307895591515771731031691422489377098686236263367916404512751010916862894129855138281961/15270994830);
5646413644023523531472659729552574911763686780871700375627426663366507837122353997075891736706811337698248660838754243486778190331522785903661556458651175061469825204821206611273990250663647381261360659950519735730925202117606150672170127523599/6;
-(8717064809960074651332043679796544474420053189621803209941154764877242575579617540185166306094852502223738126111951612525627910517501081576202264770178546608710937474005742969950212404955732473198451623024108934373488641161751158901712323446768306053/7590);
13368053158552172665210852539359893340369870651951497976111882891296650008003955172160792457229376320993686817755409436399268291095350295968657381088168219133490277914269064723832062615431730061224649980566693258603099340996988542301914519271322675688591/9366;
-(3018240015081392087620978688311925380399983229633120268872695911807562982111154053235820050168829922189401964755775948260724401542319799566237745986245598102255191922935742610508280966462644022540839619861091049093129359799053781543195492373882916779852781709/1671270);
3168221108903401670436878558215734893322849540781208738796672473984272484396317849596978630272031342024194689871467916186898192771267438982228710525079886956295106197431401217357893460897221381410667385636049264583380749631776691121592016493432807733153743581/1362;
-(1906502909997888166123201923177385833567729039151413143876271870599057704445939156075718972624235764657102074902610737729027517674632609562387841658709266014329005407533521950744449109807215808770201247724932231495252981632908042371307965561986133763291349835081839/625170);
12620737044909818561832856090355555624016887319605261762997689571062646649745107532482632213152948299491122977690702642386377706799989565320538433072623252159464115918057294711396625436506736041542403370025258477205101808063594056759977272469883621048184279331863155/3102;
-(2277782962749544647786193093844710000584561145527861043594866852071596440076587235747852022645695523676561694856919437156369375152041285055935622032497285897151359345040171955786899011718891877065760345722163063921177723102631428767500963315657655779698470644544064472947/412410);
480867224771001697116513683699011649496855159878152931805742068517626950204279499281932407966997815888727039144001177194783700618900740782637516562256421883686495287564849123342843286211825800198337962453448529082007644439295666002669973893196613894216505936316966183107269/63042;
-(321467133590936589398380572003196190798000628347443663674019204361034039315014370869884972981404460888272855773233080186485230316544246541168364468343631969225480324799028067015621769939718443419712110857965409179947456994970687005134042835903494172569465751671057323145801/30);
21954828084819816230623427376391154946628510524415478471365831649487269438565442138452375719697629394886161086199009071745032148355699097673730377547354631757000924688032674454208386076360699273002129330161098228962653466112032300056534380609136268347425707537340996253935261001/1434;
-(212191008579182060478563837945461872287372869333130175188325135660897759482730035448146388858205966593247205572842356279587190846758925659038643395344396334821348861518596112373605365460921821915409949436456394757252173113308288776951844012432992719925522001911419529928297005743854184672707/9538864545210);
4535049728315239205021018362829154800039522502549714840841943960634084990270225535043892135549408608572877513963384530325758104248133124392959997485849004663162061065909846598215984547677506961593292880902830325868627515047168286738527241360778218692535254144583771935549805772798793/138;
-(1480677868678810347062135814574727890490996459903153909612611791768134015908900253197632543925157559965099005581639883558125985134242978146873558628010545299879178299856929100217178891524159543673803785481607540954533057560554704283718320006046424881681934129216249889269701182688055001/30);
262668605206102118430195944953058387315319589613263079853189811674338644589939356676047394737425651989092115424635212785341601958591806782599345416064945113338168378810138043832891536358769753916174695061208903056552776976154338169779827414504134808495078925108569042616724875466536400029203/3486;
-(3507445546375253318342667741949651253516986113349672180095468758145505921133172244240023419466113277413785800736682458966212429373095894934752041434119711352215164597094886530278429206637066302031501412499166263642447542289785723225064437740655097573160922291075175522551466276032634749948001/30);
3053985414762198703102020975667535181880294373982570358384262712216099013630127381226064671248970168701164841859177635139561220373102316273710980757409295039231697360536736113380473515993580188204966808323201772081668695321316425977503713933894934399351348240918411922863059423248058823372368106375/16566;
-(89774288871630307843644520580598652632613559409413344033893019938898862113876440183943307724290509165539815764354254584162176909463971675189351699717189454959415964995262823025743004001915975128242680213612581748367743077774468104545056300803334120339537905702856877606246855867891786763741486735368783092237/303940201110);
2884759276650094186013272224259584433367802183845217969121665253431996536437413850716488026971480747022580182931940812035472752343853999483816456012957551120812115224018480617695398681102870535795278832047330709610952952566962354787854087985127629191094880078635434933391568643930279007460403102076327/6;
-(104203435297641800915507592803138769896333127491490532929531166632958450572751929026080533456531382568746630428805839898606131705503722712567379172835813497582706763790263297372586804372654090271412733654409742029405259156130663560412599913831518338400385290760913773448329773518728092595019422753511174189569/131070);
56081693586889462275051589447484617682656486262465693822181172785524218338951641908717980184709547885060218908068568000776767834360070429300833879662707069272458173745883316346369196270761641325176838738153606175963665383216626724969370509650128327157862121275587666040290047452303859942885949047756404247033/42;
-(5218507479961513801890596392421261361036935624312258325065379143295948300812040703848766095836974598734762472300638625802884257082786883956679824964010841565051175167717451747328911935282639583972372470105587187736495055501208701522099921363239317373617854217050435670713936357978555246779460902210809009009539232173/2291190)];

retVals=zeros(size(n));
selEven=mod(n,2)==0;
selInf=n/2+1>length(valTable);

retVals(selEven&~selInf)=valTable(1+n(selEven&~selInf)/2);
retVals(selEven&selInf)=Inf;
%The special case of n==1
retVals(n==1)=-1/2;

end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
