import React from 'react';
import {useNavigate} from 'react-router-dom';
import ToolDescription from './ToolDescription';



export default function HomeContent() {

    const navigate = useNavigate();
    const navigateToTools = () => {
        navigate('/tools');
    }

    return (
      <section className={'container home-content'}>
          <h2 className={'most-popular-tools-title'}>Most Popular Tools</h2>
          <div className={'container most-popular-tools-container'}>
              <div className={'container most-popular-tools'}>
                  <div className={'triangle-left'}></div>
                  <ToolDescription imagePath={'/images/favicon.ico'} toolTitle={'Genome Assembler'} textDescription={'Assemble a genome using short nucleotide reads'}/>
                  <ToolDescription imagePath={'/images/logo192.png'} toolTitle={'Motif Finder'} textDescription={'Find a common pattern in a list of reads'}/>
                  <ToolDescription imagePath={'/images/logo512.png'} toolTitle={'Sequence Alignment'} textDescription={'Align 2 or more reads using an alignment algorithm'}/>
                  <div className={'triangle-right'}></div>
              </div>
          </div>
          <div className={'container dot-container'}>
              <div className={'dot'}></div>
              <div className={'dot'}></div>
              <div className={'dot'}></div>
          </div>
          <button className={'blue-button see-all-tools-button'} onClick={navigateToTools}>See All Tools</button>
      </section>
    );
}