% %-----------------------------------------
 % For Molecular Physics `tMPH.bst',
 % Taylor & Francis Limited
 % generated with the docstrip utility.
 % Developed by SR Nova, Bangalore, India; modified by MKB December 2010, v.2
 % ----------------------------------------
 % *** Molecular Physics (tMPH) ***
 %
 % This bibliography style file is intended for texts in ENGLISH.
  %---------------------------------------------------------------------------
%  Produces entries of the type \bibitem{key}... in the bbl file
 % Usage of \cite is as follows:
 %   '\cite{key}' produces '[#]' in the text of the article
 %   '\cite[see][and references therein]{key}' produces '[see #, and references therein]'
 % where # is a number determined by the ordering in the References list.
 % The order in the References list is numerical.
 %---------------------------------------------------------------------

ENTRY
  { address
    author
    annote
    booktitle
    chapter
    chaptertitle
    department
    englishtitle
    edition
    editor
    eid
    howpublished
    institution
    distributor
    journal
    key
    location
    month
    note
    number
    numvols
    organization
    pages
    publisher
    school
    series
    title
    journalabbr
    transjournalabbr
    transvolume
    transyear
    transpages
    preprint
    producer
    type
    url
    volume
    year
  }
  {}
  { label extra.label sort.label short.list}
INTEGERS { output.state before.all mid.sentence after.sentence after.block }
FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}
STRINGS { s t}
FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { add.period$ write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { add.period$ " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}
FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}
FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}
FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}
FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}
FUNCTION {add.blank}
{  " " * before.all 'output.state :=
}

FUNCTION {add.comma}
{  "," * before.all 'output.state :=
}

FUNCTION {date.block}
{
  skip$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}
FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}
FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}
STRINGS {z}
FUNCTION {remove.dots}
{ 'z :=
  ""
  { z empty$ not }
  { z #1 #1 substring$
    z #2 global.max$ substring$ 'z :=
    duplicate$ "." = 'pop$
      { * }
    if$
  }
  while$
}
FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}
FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}
FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}
FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}
FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}
FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\emph{" swap$ * "}" * }
  if$
}
FUNCTION {tie.or.space.prefix}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$
}

FUNCTION {capitalize}
{ "u" change.case$ "t" change.case$ }

FUNCTION {space.word}
{ " " swap$ * " " * }
 % Here are the language-specific definitions for explicit words.
 % Each function has a name bbl.xxx where xxx is the English word.
 % The language selected here is ENGLISH
FUNCTION {bbl.and}
{ "and"}

FUNCTION {bbl.etal}
{ "{\em{et~al.}}" }

FUNCTION {bbl.editors}
{ "editors" }

FUNCTION {bbl.editor}
{ "editor" }

FUNCTION {bbl.edby}
{ "edited by" }

FUNCTION {bbl.edition}
{ "ed." }

FUNCTION {bbl.volume}
{ "Vol." }%

FUNCTION {bbl.of}
{ "of" }

FUNCTION {bbl.number}
{ "no." }

FUNCTION {bbl.nr}
{ "no." }

FUNCTION {bbl.in}
{ "in" }

FUNCTION {bbl.pages}
{ "pp." }

FUNCTION {bbl.page}
{ "p." }

FUNCTION {bbl.chapter}
{ "Chap." }

FUNCTION {bbl.techrep}
{ "" }%Tech. Rep.

FUNCTION {bbl.mthesis}
{ "Master's thesis" }

FUNCTION {bbl.phdthesis}
{ "Ph.D. thesis" }

FUNCTION {bbl.first}
{ "1st" }

FUNCTION {bbl.second}
{ "2nd" }

FUNCTION {bbl.third}
{ "3rd" }

FUNCTION {bbl.fourth}
{ "4th" }

FUNCTION {bbl.fifth}
{ "5th" }

FUNCTION {bbl.st}
{ "st" }

FUNCTION {bbl.nd}
{ "nd" }

FUNCTION {bbl.rd}
{ "rd" }

FUNCTION {bbl.plainin}
{ "in" }


FUNCTION {bbl.th}
{ "th" }

FUNCTION {bbl.unpub}
{ " (unpublished)" }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"Jun."}

MACRO {jul} {"Jul."}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sep."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

FUNCTION {eng.ord}
{ duplicate$ "1" swap$ *
  #-2 #1 substring$ "1" =
     { bbl.th * }
     { duplicate$ #-1 #1 substring$
       duplicate$ "1" =
         { pop$ bbl.st * }
         { duplicate$ "2" =
             { pop$ bbl.nd * }
             { "3" =
                 { bbl.rd * }
                 { bbl.th * }
               if$
             }
           if$
          }
       if$
     }
   if$
}

 %-------------------------------------------------------------------
 % Begin module:
 % \ProvidesFile{physjour.mbs}[2002/01/14 2.2 (PWD)]
MACRO {aa}{"Astron. \& Astrophys."}
MACRO {aasup}{"Astron. \& Astrophys. Suppl. Ser."}
MACRO {aj} {"Astron. J."}
MACRO {aph} {"Acta Phys."}
MACRO {advp} {"Adv. Phys."}
MACRO {ajp} {"Amer. J. Phys."}
MACRO {ajm} {"Amer. J. Math."}
MACRO {amsci} {"Amer. Sci."}
MACRO {anofd} {"Ann. Fluid Dyn."}
MACRO {am} {"Ann. Math."}
MACRO {ap} {"Ann. Phys. (NY)"}
MACRO {adp} {"Ann. Phys. (Leipzig)"}
MACRO {ao} {"Appl. Opt."}
MACRO {apl} {"Appl. Phys. Lett."}
MACRO {app} {"Astroparticle Phys."}
MACRO {apj} {"Astrophys. J."}
MACRO {apjsup} {"Astrophys. J. Suppl."}
MACRO {apss} {"Astrophys. Space Sci."}
MACRO {araa} {"Ann. Rev. Astron. Astrophys."}
MACRO {baas} {"Bull. Amer. Astron. Soc."}
MACRO {baps} {"Bull. Amer. Phys. Soc."}
MACRO {cmp} {"Comm. Math. Phys."}
MACRO {cpam} {"Commun. Pure Appl. Math."}
MACRO {cppcf} {"Comm. Plasma Phys. \& Controlled Fusion"}
MACRO {cpc} {"Comp. Phys. Comm."}
MACRO {cqg} {"Class. Quant. Grav."}
MACRO {cra} {"C. R. Acad. Sci. A"}
MACRO {fed} {"Fusion Eng. \& Design"}
MACRO {ft} {"Fusion Tech."}
MACRO {grg} {"Gen. Relativ. Gravit."}
MACRO {ieeens} {"IEEE Trans. Nucl. Sci."}
MACRO {ieeeps} {"IEEE Trans. Plasma Sci."}
MACRO {ijimw} {"Interntl. J. Infrared \& Millimeter Waves"}
MACRO {ip} {"Infrared Phys."}
MACRO {irp} {"Infrared Phys."}
MACRO {jap} {"J. Appl. Phys."}
MACRO {jasa} {"J. Acoust. Soc. America"}
MACRO {jcp} {"J. Comp. Phys."}
MACRO {jetp} {"Sov. Phys.--JETP"}
MACRO {jfe} {"J. Fusion Energy"}
MACRO {jfm} {"J. Fluid Mech."}
MACRO {jmp} {"J. Math. Phys."}
MACRO {jne} {"J. Nucl. Energy"}
MACRO {jnec} {"J. Nucl. Energy, C: Plasma Phys., Accelerators, Thermonucl. Res."}
MACRO {jnm} {"J. Nucl. Mat."}
MACRO {jpc} {"J. Phys. Chem."}
MACRO {jpp} {"J. Plasma Phys."}
MACRO {jpsj} {"J. Phys. Soc. Japan"}
MACRO {jsi} {"J. Sci. Instrum."}
MACRO {jvst} {"J. Vac. Sci. \& Tech."}
MACRO {nat} {"Nature"}
MACRO {nature} {"Nature"}
MACRO {nedf} {"Nucl. Eng. \& Design/Fusion"}
MACRO {nf} {"Nucl. Fusion"}
MACRO {nim} {"Nucl. Inst. \& Meth."}
MACRO {nimpr} {"Nucl. Inst. \& Meth. in Phys. Res."}
MACRO {np} {"Nucl. Phys."}
MACRO {npb} {"Nucl. Phys. B"}
MACRO {nt/f} {"Nucl. Tech./Fusion"}
MACRO {npbpc} {"Nucl. Phys. B (Proc. Suppl.)"}
MACRO {inc} {"Nuovo Cimento"}
MACRO {nc} {"Nuovo Cimento"}
MACRO {pf} {"Phys. Fluids"}
MACRO {pfa} {"Phys. Fluids A: Fluid Dyn."}
MACRO {pfb} {"Phys. Fluids B: Plasma Phys."}
MACRO {pl} {"Phys. Lett."}
MACRO {pla} {"Phys. Lett. A"}
MACRO {plb} {"Phys. Lett. B"}
MACRO {prep} {"Phys. Rep."}
MACRO {pnas} {"Proc. Nat. Acad. Sci. USA"}
MACRO {pp} {"Phys. Plasmas"}
MACRO {ppcf} {"Plasma Phys. \& Controlled Fusion"}
MACRO {phitrsl} {"Philos. Trans. Roy. Soc. London"}
MACRO {prl} {"Phys. Rev. Lett."}
MACRO {pr} {"Phys. Rev."}
MACRO {physrev} {"Phys. Rev."}
MACRO {pra} {"Phys. Rev. A"}
MACRO {prb} {"Phys. Rev. B"}
MACRO {prc} {"Phys. Rev. C"}
MACRO {prd} {"Phys. Rev. D"}
MACRO {pre} {"Phys. Rev. E"}
MACRO {ps} {"Phys. Scripta"}
MACRO {procrsl} {"Proc. Roy. Soc. London"}
MACRO {rmp} {"Rev. Mod. Phys."}
MACRO {rsi} {"Rev. Sci. Inst."}
MACRO {science} {"Science"}
MACRO {sciam} {"Sci. Am."}
MACRO {sam} {"Stud. Appl. Math."}
MACRO {sjpp} {"Sov. J. Plasma Phys."}
MACRO {spd} {"Sov. Phys.--Doklady"}
MACRO {sptp} {"Sov. Phys.--Tech. Phys."}
MACRO {spu} {"Sov. Phys.--Uspeki"}
MACRO {st} {"Sky and Telesc."}
 % End module: physjour.mbs
%% Copyright 1994-2007 Patrick W Daly
MACRO {acmcs} {"ACM Comput. Surv."}

MACRO {acta} {"Acta Inf."}

MACRO {cacm} {"Commun. ACM"}

MACRO {ibmjrd} {"IBM J. Res. Dev."}

MACRO {ibmsj} {"IBM Syst.~J."}

MACRO {ieeese} {"IEEE Trans. Software Eng."}

MACRO {ieeetc} {"IEEE Trans. Comput."}

MACRO {ieeetcad}
 {"IEEE Trans. Comput. Aid. Des."}

MACRO {ipl} {"Inf. Process. Lett."}

MACRO {jacm} {"J.~ACM"}

MACRO {jcss} {"J.~Comput. Syst. Sci."}

MACRO {scp} {"Sci. Comput. Program."}

MACRO {sicomp} {"SIAM J. Comput."}

MACRO {tocs} {"ACM Trans. Comput. Syst."}

MACRO {tods} {"ACM Trans. Database Syst."}

MACRO {tog} {"ACM Trans. Graphic."}

MACRO {toms} {"ACM Trans. Math. Software"}

MACRO {toois} {"ACM Trans. Office Inf. Syst."}

MACRO {toplas} {"ACM Trans. Progr. Lang. Syst."}

MACRO {tcs} {"Theor. Comput. Sci."}

FUNCTION {bibinfo.check}
{ swap$
  duplicate$ missing$
    {
      pop$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ pop$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
FUNCTION {bibinfo.warn}
{ swap$
  duplicate$ missing$
    {
      swap$ "missing " swap$ * " in " * cite$ * warning$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ "empty " swap$ * " in " * cite$ * warning$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
INTEGERS { nameptr namesleft numnames }


STRINGS  { bibinfo}

FUNCTION {format.names}
{ 'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{f{.}. }{vv~}{ll}{, jj}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          nameptr #2
          #2 + =
          numnames #100
          > and
            { "others" 't :=
              #1 'namesleft := }
            'skip$
          if$
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              numnames #2 >
              t "others" = not and
                { "" * }
                'skip$
              if$
              t "others" =
                {
                  " " * bbl.etal *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}

FUNCTION {format.editornames}
{ 'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{ff }{vv~}{ll}{, jj}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          nameptr #1
          #1 + =
          numnames #100
          > and
            { "others" 't :=
              #1 'namesleft := }
            'skip$
          if$
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              numnames #2 >
              t "others" = not and
                { "" * }
                'skip$
              if$
              t "others" =
                {
                  " " * bbl.etal *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.names.ed}
{
  format.editornames
}
FUNCTION {format.key}
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.authors}
{ author "author" format.names
}
FUNCTION {get.bbl.editor}
{ editor num.names$ #1 > 'bbl.editors 'bbl.editor if$ }

FUNCTION {format.editors}
{ editor "editor" format.names duplicate$ empty$ 'skip$
    {
      ", " *
      get.bbl.editor
   "" swap$ * "" *
      *
    }
  if$
}
FUNCTION {format.plaineditors}
{ editor "editor" format.names empty$ 'skip$
    {
      %"" *
      %"" *
      bbl.edby
          " " * swap$ *
          "" swap$ * "" *
    }
  if$
}

FUNCTION {format.url}
{  url empty$
    'skip$
    { add.period$ add.blank " $<${" url * "}$>$" * output }
  if$
 }
 
 FUNCTION {format.unpuburl}
{  url empty$
    'skip$
    { add.comma add.blank " $<${" url * "}$>$" * output }
  if$
 }

FUNCTION {format.note}
{ note empty$
    { "" }
    { note #1 #1 substring$
      duplicate$ "{" =
        'skip$
        { output.state mid.sentence =
          { "t" }
          { "t" }
        if$
        change.case$
        }
      if$
      note #2 global.max$ substring$ * "note" bibinfo.check
    }
  if$
}

FUNCTION {format.unpubnote}
{ note empty$
    {  bbl.unpub  *  }
    {  note "note" bibinfo.check
       bbl.unpub  * swap$  " "   }
      if$
    }
FUNCTION {end.quote.title}
{ title empty$
    'skip$
    { before.all 'output.state := }
  if$
}
FUNCTION {format.trnbookedition}
{ edition duplicate$ empty$ 'skip$
    { "edition" bibinfo.check
      add.blank ", " swap$ *
      "" *
    }
  if$
}
FUNCTION {format.title}
{ title
  duplicate$ empty$ 'skip$
    {  }%"t" change.case$
  if$
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      %emphasize
    }
  if$
}
FUNCTION {format.trntitle}
{ title empty$
    { "" }
    {add.blank title "title" bibinfo.check
      duplicate$ empty$ 'skip$
      { emphasize
       "" swap$ *
      "" *
      }
    if$
    }
  if$
}


FUNCTION {format.englishtitle}
{ englishtitle empty$
    { "" }
    { add.blank englishtitle "englishtitle" bibinfo.check
      duplicate$ empty$ 'skip$
        { emphasize
           "[" swap$ *
          "]" *
        }
      if$
    }
  if$
}

FUNCTION {format.chaptertitle}
{ chaptertitle empty$
    { "" }
    {chaptertitle "chaptertitle" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      ". " *
    }
  if$
}

FUNCTION {format.address}
{ address empty$
    { "" }
    {address "address" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      ", " *
    }
  if$
}

FUNCTION {format.producer}
{ producer empty$
    { "" }
    {producer "producer" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { "film produced by " swap$ *
      ", " *
    }
  if$
}

FUNCTION {format.distributor}
{ distributor empty$
    { "" }
    {distributor "distributor" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { "(distributed by " swap$ *
      "; " *
    }
  if$
}

FUNCTION {format.institution}
{ institution empty$
    { "" }
    {institution "institution" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      "" *
    }
  if$
}

FUNCTION {format.thesisinstitution}
{ institution empty$
    { "" }
    {institution "institution" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      "," *
    }
  if$
}

FUNCTION {format.unpubinstitution}
{ institution empty$
    { "" }
    {institution "institution" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      "," *
    }
  if$
}

FUNCTION {format.department}
{ department empty$
    { " " }
    {department "department" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      "" *
    }
  if$
}

FUNCTION {format.school}
{ school empty$
    { "" }
    {school "school" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      "" *
    }
  if$
}

FUNCTION {format.transjournalabbr}
{ transjournalabbr empty$
    { "" }
    {transjournalabbr "transjournalabbr" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " [" swap$ *
      "" *
    }
  if$
}

FUNCTION {format.transvolume}
{ transvolume empty$
    { "" }
    {transvolume "transvolume" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      "" *
    }
  if$
}

FUNCTION {format.transyear}
{ transyear empty$
    { "" }
    {transyear "transyear" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { "(" swap$ *
      ")]" *
    }
  if$
}

FUNCTION {format.transpages}
{ transpages empty$
    { "" }
    {transpages "transpages" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      " " *
    }
  if$
}
FUNCTION {format.journaltitle}
{ journalabbr empty$
    { journal "journal" bibinfo.check }
    {journalabbr "journalabbr" bibinfo.check}
  if$
  duplicate$ empty$
    'skip$
    { " " swap$ *
      " " *
    }
  if$
}

FUNCTION {format.plaintitle}
{ title
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "" swap$ *
      ". " *
    }
  if$
}
FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {n.dashify}
{
  't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {unen.dashify}
{
  't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "" *
              t "" global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "" *
                  t "" global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {word.in}
{ bbl.in
  " " * }

FUNCTION {word.plainin}
{ bbl.plainin
  " " * }
  
FUNCTION {word.unpublished}
{ bbl.unpub
  " " * }


FUNCTION {format.date}
{
  ""
  duplicate$ empty$
  year  "year"  bibinfo.check duplicate$ empty$
    { swap$ 'skip$
        { "there's a month but no year in " cite$ * warning$ }
      if$
      *
    }
    { swap$ 'skip$
        {
          swap$
          " " * swap$
        }
      if$
      *
      remove.dots
    }
  if$
  duplicate$ empty$
    'skip$
    {
      before.all 'output.state :=
    " " swap$ * "" *
    }
  if$
}

FUNCTION {format.unpubdate}
{
  ""
  duplicate$ empty$
  year  "year"  bibinfo.check duplicate$ empty$
    { swap$ 'skip$
        { "there's a month but no year in " cite$ * warning$ }
      if$
      *
    }
    { swap$ 'skip$
        {
          swap$
          " " * swap$
        }
      if$
      *
      remove.dots
    }
  if$
  duplicate$ empty$
    'skip$
    {
      before.all 'output.state :=
    " " swap$ * "" *
    }
  if$
}

FUNCTION {format.fdate}
{ year "year" bibinfo.check duplicate$ empty$
    {
      "empty year in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  " released " swap$ * ")" *
}


FUNCTION {format.dbasedate}
{
  ""
  duplicate$ empty$
  year  "year"  bibinfo.check duplicate$ empty$
    { swap$ 'skip$
        { "there's a month but no year in " cite$ * warning$ }
      if$
      *
    }
    { swap$ 'skip$
        {
          swap$
          " " * swap$
        }
      if$
      *
      remove.dots
    }
  if$
  duplicate$ empty$
    'skip$
    {
      before.all 'output.state :=
    ", " swap$ * "" *
    }
  if$
}


FUNCTION {format.plaindate}
{ year "year" bibinfo.check duplicate$ empty$
    {
      "empty year in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  " (" swap$ * ")" *
}
FUNCTION{format.year}
{ year "year" bibinfo.check duplicate$ empty$
    {  "empty year in " cite$ *
       warning$
    }
    {
  " (" swap$ * ")" *
    }
  if$
}
FUNCTION {format.btitle}
{ title "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      emphasize
    }
  if$
}


FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}
FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.prefix
      "volume" bibinfo.check * *
      series "series" bibinfo.check
      duplicate$ empty$ 'pop$
        { emphasize ", " * swap$ * }
      if$
      "volume and number" number either.or.check
    }
  if$
}
FUNCTION {format.number.series}
{ volume empty$
    { number empty$
        { series field.or.null }
        { series empty$
            { number "number" bibinfo.check }
            { output.state mid.sentence =
                { bbl.number }
                { bbl.number capitalize }
              if$
              number tie.or.space.prefix "number" bibinfo.check * *
              bbl.in space.word *
              series "series" bibinfo.check *
            }
          if$
        }
      if$
    }
    { "" }
  if$
}
FUNCTION {is.num}
{ chr.to.int$
  duplicate$ "0" chr.to.int$ < not
  swap$ "9" chr.to.int$ > not and
}

FUNCTION {extract.num}
{ duplicate$ 't :=
  "" 's :=
  { t empty$ not }
  { t #1 #1 substring$
    t #2 global.max$ substring$ 't :=
    duplicate$ is.num
      { s swap$ * 's := }
      { pop$ "" 't := }
    if$
  }
  while$
  s empty$
    'skip$
    { pop$ s }
  if$
}

FUNCTION {convert.edition}
{ extract.num "l" change.case$ 's :=
  s "first" = s "1" = or
    { bbl.first 't := }
    { s "second" = s "2" = or
        { bbl.second 't := }
        { s "third" = s "3" = or
            { bbl.third 't := }
            { s "fourth" = s "4" = or
                { bbl.fourth 't := }
                { s "fifth" = s "5" = or
                    { bbl.fifth 't := }
                    { s #1 #1 substring$ is.num
                        { s eng.ord 't := }
                        { edition 't := }
                      if$
                    }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
  t
}
FUNCTION {format.number.plainseries}
{ volume empty$
    { number empty$
        { series field.or.null }
        { series empty$
            { number "number" bibinfo.check }
            { output.state mid.sentence =
                { bbl.number }
                { bbl.number capitalize }
              if$
              number tie.or.space.prefix "number" bibinfo.check * *
              bbl.plainin space.word *
              series "series" bibinfo.check *
            }
          if$
        }
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition duplicate$ empty$ 'skip$
    {
      convert.edition
      output.state mid.sentence =
        { "l" }
        { "t" }
      if$ change.case$
      "edition" bibinfo.check
      " " * bbl.edition *
    }
  if$
}
INTEGERS { multiresult }
FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}
FUNCTION {format.pages}
{ pages duplicate$ empty$ 'skip$
    { duplicate$ multi.page.check
        {
          bbl.pages swap$
          n.dashify
        }
        {
          bbl.page swap$
        }
      if$
      tie.or.space.prefix
      "pages" bibinfo.check
      * *
    }
  if$
}


FUNCTION {format.journal.pages}
{ pages empty$ 'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        {
          ", " *
          pages n.dashify *
        }
      if$
    }
  if$
}

FUNCTION {format.journal.artpages}
{ pages empty$ 'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        {
          ", " *
          pages unen.dashify *
        }
      if$
    }
  if$
}


FUNCTION {format.journal.eid}
{ eid "eid" bibinfo.check
  duplicate$ empty$ 'pop$
    { swap$ duplicate$ empty$ 'skip$
      {
          ", " *
      }
      if$
      swap$ *
    }
  if$
}
FUNCTION {bold}
{ duplicate$ empty$
    { pop$ "" }
    { "\textbf{" swap$ * "}" * }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
      bold "volume" bibinfo.check
    }
  if$
  number "number" bibinfo.check duplicate$ empty$ 'skip$
    {
      swap$ duplicate$ empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
      swap$
      " (" swap$ * ")" *
    }
  if$ *
  eid empty$
    { format.journal.pages }
    { format.journal.eid }
  if$
}

FUNCTION {format.vol.num.artpages}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
      bold "volume" bibinfo.check
    }
  if$
  number "number" bibinfo.check duplicate$ empty$ 'skip$
    {
      swap$ duplicate$ empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
      swap$
      " (" swap$ * ")" *
    }
  if$ *
  eid empty$
    { format.journal.artpages }
    { format.journal.eid }
  if$
}


FUNCTION {format.volume}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
     bold  "volume" bibinfo.check
      "" *
    }
  if$
}

FUNCTION {format.numvols}
{ numvols field.or.null
  duplicate$ empty$ 'skip$
    {
     "numvols" bibinfo.check
      " vols." *
    }
  if$
}
FUNCTION {format.chapter.pages}
{ chapter empty$
    { "" }
    { type empty$
        { bbl.chapter }
        { type "l" change.case$
          "type" bibinfo.check
        }
      if$
      chapter tie.or.space.prefix
      "chapter" bibinfo.check
      * *
    }
  if$
}

FUNCTION {format.booktitle}
{
  booktitle "booktitle" bibinfo.check
  emphasize
}
FUNCTION {format.in.booktitle}
{
  booktitle "booktitle" bibinfo.check
  emphasize
      word.in swap$ *
}
FUNCTION {format.in.ed.booktitle}
{ format.booktitle duplicate$ empty$ 'skip$
    {
      editor "editor"  format.names.ed duplicate$ empty$ 'pop$
        {
          swap$  "," *
          " " * swap$ *
          swap$
          "," *
          " " * swap$
          * }
      if$
      word.in swap$ *
    }
  if$
}
FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}
FUNCTION {format.thesis.type}
{ type duplicate$ empty$
    'pop$
    { swap$ pop$
      %"t" change.case$ "type" bibinfo.check
    }
  if$
}
FUNCTION {format.tr.number}
{ number "number" bibinfo.check
  type duplicate$ empty$
    { pop$ bbl.techrep }
    'skip$
  if$
  "type" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
}
FUNCTION {format.article.crossref}
{
  key duplicate$ empty$
    { pop$
      journal duplicate$ empty$
        { "need key or journal for " cite$ * " to crossref " * crossref * warning$ }
        { "journal" bibinfo.check emphasize word.in swap$ * }
      if$
    }
    { word.in swap$ * " " *}
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  "editor" bibinfo.check
  editor num.names$ duplicate$
  #2 >
    { pop$
      "editor" bibinfo.check
      " " * bbl.etal
      *
    }
    { #2 <
        'skip$
        { editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            {
              "editor" bibinfo.check
              " " * bbl.etal
              *
            }
            {
             bbl.and space.word
              * editor #2 "{vv~}{ll}" format.name$
              "editor" bibinfo.check
              *
            }
          if$
        }
      if$
    }
  if$
}
FUNCTION {format.book.crossref}
{ volume duplicate$ empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      pop$ word.in
    }
    { bbl.volume
      swap$ tie.or.space.prefix "volume" bibinfo.check * * bbl.of space.word *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { series empty$
            { "need editor, key, or series for " cite$ * " to crossref " *
              crossref * warning$
              "" *
            }
            { series emphasize * }
          if$
        }
        { key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.incoll.inproc.crossref}
{
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { format.booktitle duplicate$ empty$
            { "need editor, key, or booktitle for " cite$ * " to crossref " *
              crossref * warning$
            }
            { word.in swap$ * }
          if$
        }
        { word.in key * " " *}
      if$
    }
    { word.in format.crossref.editor * " " *}
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.org.or.pub}
{ 't :=
  " ("
  address empty$ t empty$ and
    'skip$
    {
      t empty$
        'skip$
        { address empty$
            'skip$
            { t * }
          if$
          "" swap$ *  ", " *
        }
      if$
      address "address" bibinfo.check *
    }
  if$
   * "" swap$
  year "year" output.check
  note empty$
  { "" }
  {
  add.blank
  swap$ * ", "  swap$ * "\unskip; "
  note duplicate$ empty$
    { "" }
    {* "" }
  if$
}
  if$
  ")"  *  *
}


FUNCTION {noformat.org.or.pub}
{ 't :=
  ""
  address empty$ t empty$ and
    'skip$
    {
      address "address" bibinfo.check *
      t empty$
        'skip$
        { address empty$
            'skip$
            { ": " * }
          if$
          t *
        }
      if$
    }
  if$
}

FUNCTION {format.publisher.address}
{ publisher empty$
    { publisher "publisher" bibinfo.warn noformat.org.or.pub }
    { add.blank
      publisher "publisher" bibinfo.warn format.org.or.pub }
  if$
}

FUNCTION {format.organization.address}
{ organization "organization" bibinfo.check format.org.or.pub
}

FUNCTION {format.location}
{ location empty$
    { "" }
    { location "location" bibinfo.check}
  if$
}
FUNCTION {format.month}
{ month empty$
    { "" }
    { month "month" bibinfo.check}
  if$
}

FUNCTION {format.inproc.booktitle}
{
  booktitle "booktitle" bibinfo.check
       " " word.plainin swap$ emphasize * swap$ emphasize *%Proceedings of the
}

FUNCTION {format.plainin.booktitle}
{
  booktitle "booktitle" bibinfo.check
}

STRINGS {oldname}

FUNCTION {name.or.dash}
{ 's :=
   oldname empty$
     { s 'oldname := s }
     { s oldname =
         { "---------" }
         { s 'oldname := s }
       if$
     }
   if$
}

FUNCTION {article}
{ output.bibitem
  format.authors output
  crossref missing$
    {
      format.journaltitle output
      end.quote.title
      add.blank
      format.vol.num.artpages output
      format.plaindate "year" output.check
    }
    { format.article.crossref output.nonnull
      format.journal.artpages output
    }
  if$
  format.note output
  fin.entry
}


FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      format.numvols output
      format.edition output
      end.quote.title
      add.blank
      format.number.series output
      format.publisher.address output
    }
    {
      format.book.crossref output.nonnull
    }
  if$
  format.pages output
  format.note output
  fin.entry
}
FUNCTION {booklet}
{ output.bibitem
  format.authors output
  format.title "title" output.check
  howpublished "howpublished" bibinfo.check output
  address "address" bibinfo.check output
  format.date output
  format.note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  format.authors "author" output.check
  crossref missing$
  {format.in.booktitle "booktitle" output.check
  editor "editor" format.names.ed
  format.plaineditors "author and editor" output.check
   format.publisher.address output
      format.bvolume output
      format.chapter.pages "chapter and pages" output.check
      format.edition output
      format.number.series output
    }
    {
      format.chapter.pages "chapter and pages" output.check
      format.book.crossref output.nonnull
    }
  if$
  %format.date "year" output.check
  format.pages "pages" output.check
  format.note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  crossref missing$
    {
      format.in.booktitle "booktitle" output.check
      editor "editor" format.names.ed
      format.plaineditors "author and editor" output.check
      format.bvolume output
      format.edition output
      format.chapter.pages output
      format.number.series output
      format.publisher.address output
      %format.date "year" output.check
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  format.pages output
  format.note output
  fin.entry
}
FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  %format.plaindate "year" output.check
  format.plaintitle "title" output.check
  end.quote.title
  crossref missing$
    { format.inproc.booktitle "booktitle" output.check
      format.edition output
      editor empty$
        { }
        { editor "editor" format.names.ed
          format.plaineditors "author and editor" output.check
        }
      if$
      format.bvolume output
      format.number.plainseries output
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.location output
  format.month output
  format.note output
  format.publisher.address output
  format.chapter.pages output
  format.pages "pages" output.check
  fin.entry
}
FUNCTION {conference} { inproceedings }
FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization "organization" bibinfo.check
      duplicate$ empty$ 'pop$
        { output
          address "address" bibinfo.check output
        }
      if$
    }
    { format.authors output.nonnull }
  if$
  format.btitle "title" output.check
  author empty$
    { organization empty$
        {
          address "address" bibinfo.check output
        }
        'skip$
      if$
    }
    {
      organization "organization" bibinfo.check output
      address "address" bibinfo.check output
    }
  if$
  format.edition output
  format.date output
  format.note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  bbl.mthesis format.thesis.type output.nonnull
  school "school" bibinfo.warn output
  format.thesisinstitution "institution" output.check
  address "address" bibinfo.check output
  format.date "year" output.check
  format.note output
  fin.entry
}

FUNCTION {database}
{ output.bibitem
  format.authors output
  format.title output
  format.institution output
  format.department output
  format.note output
  format.dbasedate output
  end.quote.title
  format.url "url" bibinfo.check output
  fin.entry
}


FUNCTION {misc}
{ output.bibitem
  format.authors output
  format.title output
  end.quote.title
  add.blank
  %format.publisher.address output
  publisher "publisher" bibinfo.check output
  format.address output
  howpublished "howpublished" bibinfo.check output
  format.date output
  format.note output
  format.month output
  fin.entry
  empty.misc.check
}

FUNCTION {USPATENT}
{ output.bibitem
  format.authors output
  format.title output
  format.note output
  end.quote.title
  add.blank
  publisher "publisher" bibinfo.check output
  address "address" bibinfo.check output
  format.month output
   format.date output
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  bbl.phdthesis format.thesis.type output.nonnull
  format.department "department" bibinfo.check output
  format.school "school" bibinfo.check output
  format.thesisinstitution "institution" output.check
  end.quote.title
  format.address "address" bibinfo.check output
  end.quote.title
  format.date "year" output.check
  format.note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization "organization" bibinfo.check output
    }
    { format.editors output.nonnull }
  if$
  format.btitle "title" output.check
  format.bvolume output
  editor empty$
    { publisher empty$
        {  format.number.series output }
        {
          format.number.series output
          format.publisher.address output
        }
      if$
    }
    { publisher empty$
        {
          format.number.series output
          format.organization.address output }
        {
          format.number.series output
          organization "organization" bibinfo.check output
          format.publisher.address output
        }
      if$
     }
  if$
      format.date "year" output.check
  format.note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  format.tr.number output.nonnull
  institution "institution" bibinfo.check output
  address "address" bibinfo.check output
   format.date output
  end.quote.title
  add.blank
  format.note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  format.unpubinstitution "institution" output.check
  end.quote.title
  format.address "address" output.check
  format.unpubdate output
  end.quote.title
  format.unpubnote output
  new.block
  format.unpuburl "url" bibinfo.check output
  fin.entry
}

FUNCTION {EPRINT}
{ output.bibitem
  format.authors "author" output.check
  format.title output
  preprint "preprint" bibinfo.check output
  format.plaindate "year" output.check
  format.note output
  format.url "url" bibinfo.check output
  fin.entry
}

FUNCTION {film}
{ output.bibitem
  format.authors "author" output.check
  format.btitle "title" output.check
  format.producer "producer" output.check
  end.quote.title
  address "address" bibinfo.check output
  end.quote.title
  add.blank
  format.distributor "distributor" output.check
  format.fdate "year" output.check
  format.note "note" output.check
  fin.entry
}


FUNCTION {translatedbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
      add.blank
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
      }
  if$
  format.plaindate "year" output.check
  format.chaptertitle "chaptertitle" output.check
  format.trntitle "title" output.check
  format.englishtitle "englishtitle" output.check
  format.trnbookedition "edition" output.check
  crossref missing$
    { format.bvolume output
      format.number.series output
      format.publisher.address output
    }
    {
      format.book.crossref output.nonnull
    }
  if$
  format.note output
  fin.entry
}


FUNCTION {translatedarticle}
{ output.bibitem
  format.authors "author" output.check
  publisher "publisher" bibinfo.check output
  journalabbr "journalabbr" output.check
  end.quote.title
  add.blank
  format.vol.num.pages output
  %format.volume output
  %format.pages output
  format.plaindate "year" output.check
  end.quote.title
  format.transjournalabbr "transjournalabbr" output.check
  end.quote.title
  format.transvolume "transvolume" output.check
  format.transpages "transpages" output.check
  end.quote.title
  format.transyear "transyear" output.check
  format.note output
  fin.entry
}

FUNCTION {default.type} { misc }
READ
FUNCTION {sortify}
{ purify$
  "l" change.case$
}
INTEGERS { len }
FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}
FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}"
      format.name$ 't :=
      nameptr #1 >
        {
          nameptr #1
          #1 + =
          numnames #3
          > and
            { "others" 't :=
              #1 'namesleft := }
            'skip$
          if$
          "   "  *
          namesleft #1 = t "others" = and
            { "zzzzz" * }
            { t sortify * }
          if$
        }
        { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}
FUNCTION {author.sort}
{ author empty$
    { key empty$
        { "to sort, need author or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { key empty$
            { "to sort, need author, editor, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
        { key empty$
            { "to sort, need author, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
        { key empty$
            { "to sort, need editor, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}
FUNCTION {presort}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
        'editor.organization.sort
        { type$ "manual" =
            'author.organization.sort
            'author.sort
          if$
        }
      if$
    }
  if$
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}
ITERATE {presort}
%SORT
STRINGS { longest.label }
INTEGERS { number.label longest.label.width }
FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}
FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}
%EXECUTE {initialize.longest.label}
ITERATE {longest.label.pass}
FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{" number.label int.to.str$ * "}" *
  write$ newline$
  "\providecommand{\url}[1]{\texttt{#1}}"
  write$ newline$
  "\providecommand{\urlprefix}{URL }"
  write$ newline$
  "\markboth{Taylor \& Francis and I.T. Consultant}{Molecular Physics}"
  write$ newline$
}
EXECUTE {begin.bib}
EXECUTE {init.state.consts}
ITERATE {call.type$}
FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}
EXECUTE {end.bib}
%% End of customized bst file
%%
%% End of file `tMPH2.bst'.
