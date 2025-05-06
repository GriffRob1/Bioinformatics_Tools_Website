import React from 'react';
import TitleArea from "../components/TitleArea";
import HomeContent from "../components/HomeContent";

export default function Home({toolsList}) {
    return (
        <div className='home-page'>
            <TitleArea />
            <HomeContent toolsList={toolsList}/>
        </div>
    );
}