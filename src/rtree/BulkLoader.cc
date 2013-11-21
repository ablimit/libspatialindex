/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

#include <cstring>
#include <stdio.h>
#include <cmath>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <spatialindex/SpatialIndex.h>

#include "RTree.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"

using namespace SpatialIndex;
using namespace SpatialIndex::RTree;

//
// ExternalSorter::Record
//
    ExternalSorter::Record::Record()
: m_pData(0)
{
}

    ExternalSorter::Record::Record(const Region& r, id_type id, uint32_t len, byte* pData, uint32_t s)
: m_r(r), m_id(id), m_len(len), m_pData(pData), m_s(s)
{
}

ExternalSorter::Record::~Record()
{
    delete[] m_pData;
}

bool ExternalSorter::Record::operator<(const Record& r) const
{
    if (m_s != r.m_s)
	throw Tools::IllegalStateException("ExternalSorter::Record::operator<: Incompatible sorting dimensions.");

    if (m_r.m_pHigh[m_s] + m_r.m_pLow[m_s] < r.m_r.m_pHigh[m_s] + r.m_r.m_pLow[m_s])
	return true;
    else
	return false;
}

void ExternalSorter::Record::storeToFile(Tools::TemporaryFile& f)
{
    f.write(static_cast<uint64_t>(m_id));
    f.write(m_r.m_dimension);
    f.write(m_s);

    for (uint32_t i = 0; i < m_r.m_dimension; ++i)
    {
	f.write(m_r.m_pLow[i]);
	f.write(m_r.m_pHigh[i]);
    }

    f.write(m_len);
    if (m_len > 0) f.write(m_len, m_pData);
}

void ExternalSorter::Record::loadFromFile(Tools::TemporaryFile& f)
{
    m_id = static_cast<id_type>(f.readUInt64());
    uint32_t dim = f.readUInt32();
    m_s = f.readUInt32();

    if (dim != m_r.m_dimension)
    {
	delete[] m_r.m_pLow;
	delete[] m_r.m_pHigh;
	m_r.m_dimension = dim;
	m_r.m_pLow = new double[dim];
	m_r.m_pHigh = new double[dim];
    }

    for (uint32_t i = 0; i < m_r.m_dimension; ++i)
    {
	m_r.m_pLow[i] = f.readDouble();
	m_r.m_pHigh[i] = f.readDouble();
    }

    m_len = f.readUInt32();
    delete[] m_pData; m_pData = 0;
    if (m_len > 0) f.readBytes(m_len, &m_pData);
}

//
// ExternalSorter
//
    ExternalSorter::ExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages)
: m_bInsertionPhase(true), m_u32PageSize(u32PageSize),
    m_u32BufferPages(u32BufferPages), m_u64TotalEntries(0), m_stI(0), last_dim(-1)
{
}

ExternalSorter::~ExternalSorter()
{
    for (m_stI = 0; m_stI < m_buffer.size(); ++m_stI) delete m_buffer[m_stI];
}

void ExternalSorter::insert(Record* r)
{
    if (m_bInsertionPhase == false)
	throw Tools::IllegalStateException("ExternalSorter::insert: Input has already been sorted.");

    m_buffer.push_back(r);
    ++m_u64TotalEntries;

    // this will create the initial, sorted buckets before the
    // external merge sort.
    //if (m_buffer.size() >= static_cast<uint64_t>m_u32PageSize * static_cast<uint64_t>m_u32BufferPages)
    if (m_buffer.size() >= 100000000) // 100 Million objects 
    {
	std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
	Tools::TemporaryFile* tf = new Tools::TemporaryFile();
	for (size_t j = 0; j < m_buffer.size(); ++j)
	{
	    m_buffer[j]->storeToFile(*tf);
	    delete m_buffer[j];
	}
	m_buffer.clear();
	tf->rewindForReading();
	m_runs.push_back(Tools::SmartPointer<Tools::TemporaryFile>(tf));
    }
}

void ExternalSorter::split(uint32_t K,uint32_t dim, Region &r, std::vector<Record*> &node)
{
    double c1 [] = {0.0, 0.0};
    double c2 [] = {0.0, 0.0};
    uint32_t adim = (dim+1)%2 ; // another dimension 
    Region p(c1,c2,2);

    std::vector<Record*>::size_type len = 0;

    if (getTotalEntries() >K)
    {
	if (dim != last_dim)
	    sort(dim,K);
	len = K;
	p = m_buffer[K-1]->m_r ; 

	c1[dim] = p.getHigh(dim);
	c1[adim] = universe.getLow(adim);
	c2[dim] = c1[dim];
	c2[adim] = universe.getHigh(adim);

	memcpy(r.m_pLow, universe.m_pLow,   2 * sizeof(double));
	memcpy(r.m_pHigh, c2, 2 * sizeof(double));

	memcpy(universe.m_pLow, c1, 2 * sizeof(double));

    }
    else if (getTotalEntries() <= K)
    {
	//for (len = 0 ; len < m_buffer.size(); len++) // temporarily used len as a index var
	//    r.combineRegion(m_buffer[len]->m_r);
	r = universe; 
	len = getTotalEntries();
    }

    // update container 
    std::vector<Record*>::iterator it=m_buffer.begin(); 
    //node.assign(it,it+len);
    m_buffer.erase(it,it+len);
    m_u64TotalEntries = m_buffer.size();
}

float ExternalSorter::getCost(uint32_t K, uint32_t dim)
{
    uint32_t adim = (dim+1)%2 ; // another dimension 
    float cost = 0.0; 
    //TODO:
    //this code should be optimzed so that we only need to sort top alpha*K elements ; 
    sort(dim);
    Region kr = m_buffer[K-1]->m_r;
    double c1 [] = {0.0, 0.0};
    double c2 [] = {0.0, 0.0};

    c1[dim] = kr.getHigh(dim);
    c1[adim] = universe.getLow(adim);
    c2[dim] = c1[dim];
    c2[adim] = universe.getHigh(adim);
    LineSegment lseg(c1,c2,universe.getDimension());

    // iterate left side
    for (std::vector<Record*>::size_type i = K; i >0 ;i--)
    {
	if (m_buffer[i]->m_r.intersectsLineSegment(lseg))
	    cost += 1.0 ;
	else 
	{
	    if (i == 0 && m_buffer[i]->m_r.intersectsLineSegment(lseg))
		cost += 1.0 ;
	    else
		break;
	}
    }
    // iterate right side 
    for (std::vector<Record*>::size_type i = K+1; i <getTotalEntries(); i++)
    {
	if (m_buffer[i]->m_r.intersectsLineSegment(lseg))
	    cost += 1.0 ;
	else 
	    break;
    }
    return cost;
}

Region ExternalSorter::getUniverse(){
    //init
    if (getTotalEntries()>0)
	universe = m_buffer[0]->m_r;

    for (std::vector<Record*>::iterator it=m_buffer.begin(); it != m_buffer.end(); ++it)
    {
	universe.combineRegion((*it)->m_r);
    }
    return universe;
}

void ExternalSorter::sort(uint32_t dim, uint64_t K)
{
    if (K > 0 ) {
	switch(dim){
	    case 0 :
		std::partial_sort(m_buffer.begin(), m_buffer.begin()+K, m_buffer.end(), Record::SortAscendingX());
		break;

	    case 1:
		std::partial_sort(m_buffer.begin(), m_buffer.begin()+K, m_buffer.end(), Record::SortAscendingY());
		break;
	    default:
		throw Tools::IllegalStateException("ExternalSorter::sort Incompatible sorting dimensions.");
	}
    }
    else {
	switch(dim){
	    case 0 :
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscendingX());
		break;

	    case 1:
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscendingY());
		break;

	    default:
		throw Tools::IllegalStateException("ExternalSorter::sort Incompatible sorting dimensions.");
	}
    }

    last_dim = dim ;
}

void ExternalSorter::sort()
{
    if (m_bInsertionPhase == false)
	throw Tools::IllegalStateException("ExternalSorter::sort: Input has already been sorted.");

    if (m_runs.empty())
    {
	// The data fits in main memory. No need to store to disk.
	std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
	m_bInsertionPhase = false;
	return;
    }

    if (m_buffer.size() > 0)
    {
	// Whatever remained in the buffer (if not filled) needs to be stored
	// as the final bucket.
	std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
	Tools::TemporaryFile* tf = new Tools::TemporaryFile();
	for (size_t j = 0; j < m_buffer.size(); ++j)
	{
	    m_buffer[j]->storeToFile(*tf);
	    delete m_buffer[j];
	}
	m_buffer.clear();
	tf->rewindForReading();
	m_runs.push_back(Tools::SmartPointer<Tools::TemporaryFile>(tf));
    }

    if (m_runs.size() == 1)
    {
	m_sortedFile = m_runs.front();
    }
    else
    {
	Record* r = 0;

	while (m_runs.size() > 1)
	{
	    Tools::SmartPointer<Tools::TemporaryFile> tf(new Tools::TemporaryFile());
	    std::vector<Tools::SmartPointer<Tools::TemporaryFile> > buckets;
	    std::vector<std::queue<Record*> > buffers;
	    std::priority_queue<PQEntry, std::vector<PQEntry>, PQEntry::SortAscending> pq;

	    // initialize buffers and priority queue.
	    std::list<Tools::SmartPointer<Tools::TemporaryFile> >::iterator it = m_runs.begin();
	    for (uint32_t i = 0; i < (std::min)(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages); ++i)
	    {
		buckets.push_back(*it);
		buffers.push_back(std::queue<Record*>());

		r = new Record();
		r->loadFromFile(**it);
		// a run cannot be empty initially, so this should never fail.
		pq.push(PQEntry(r, i));

		for (uint32_t j = 0; j < m_u32PageSize - 1; ++j)
		{
		    // fill the buffer with the rest of the page of records.
		    try
		    {
			r = new Record();
			r->loadFromFile(**it);
			buffers.back().push(r);
		    }
		    catch (Tools::EndOfStreamException)
		    {
			delete r;
			break;
		    }
		}
		++it;
	    }

	    // exhaust buckets, buffers, and priority queue.
	    while (! pq.empty())
	    {
		PQEntry e = pq.top(); pq.pop();
		e.m_r->storeToFile(*tf);
		delete e.m_r;

		if (! buckets[e.m_u32Index]->eof() && buffers[e.m_u32Index].empty())
		{
		    for (uint32_t j = 0; j < m_u32PageSize; ++j)
		    {
			try
			{
			    r = new Record();
			    r->loadFromFile(*buckets[e.m_u32Index]);
			    buffers[e.m_u32Index].push(r);
			}
			catch (Tools::EndOfStreamException)
			{
			    delete r;
			    break;
			}
		    }
		}

		if (! buffers[e.m_u32Index].empty())
		{
		    e.m_r = buffers[e.m_u32Index].front();
		    buffers[e.m_u32Index].pop();
		    pq.push(e);
		}
	    }

	    tf->rewindForReading();

	    // check if another pass is needed.
	    uint32_t u32Count = std::min(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages);
	    for (uint32_t i = 0; i < u32Count; ++i)
	    {
		m_runs.pop_front();
	    }

	    if (m_runs.size() == 0)
	    {
		m_sortedFile = tf;
		break;
	    }
	    else
	    {
		m_runs.push_back(tf);
	    }
	}
    }

    m_bInsertionPhase = false;
}

ExternalSorter::Record* ExternalSorter::getNextRecord()
{
    if (m_bInsertionPhase == true)
	throw Tools::IllegalStateException("ExternalSorter::getNextRecord: Input has not been sorted yet.");

    Record* ret;

    if (m_sortedFile.get() == 0)
    {
	if (m_stI < m_buffer.size())
	{
	    ret = m_buffer[m_stI];
	    m_buffer[m_stI] = 0;
	    ++m_stI;
	}
	else
	    throw Tools::EndOfStreamException("");
    }
    else
    {
	ret = new Record();
	ret->loadFromFile(*m_sortedFile);
    }

    return ret;
}

inline uint64_t ExternalSorter::getTotalEntries() const
{
    return m_u64TotalEntries;
}

// BulkLoader
void BulkLoader::bulkLoadUsingSTRIP(
	SpatialIndex::RTree::RTree* pTree,
	IDataStream& stream,
	uint32_t b,
	uint32_t dim
	) {
    if (! stream.hasNext())
	throw Tools::IllegalArgumentException(
		"RTree::BulkLoader::bulkLoadUsingRplus: Empty data stream given."
		);
    NodePtr n = pTree->readNode(pTree->m_rootID);
    pTree->deleteNode(n.get());

#ifndef NDEBUG
    std::cerr << "RTree::BulkLoader:STRIP Sorting data." << std::endl;
#endif

    Tools::SmartPointer<ExternalSorter> es = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(10000,10000));

    while (stream.hasNext())
    {
	Data* d = reinterpret_cast<Data*>(stream.getNext());
	if (d == 0)
	    throw Tools::IllegalArgumentException(
		    "bulkLoadUsingSTRIP: RTree bulk load expects SpatialIndex::RTree::Data entries."
		    );

	es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, dim));
	d->m_pData = 0;
	delete d;
    }

    // sort by dim 
    es->sort(dim);

    pTree->m_stats.m_u64Data = es->getTotalEntries();

    // create index levels.
    uint32_t level = 0;
    std::vector<ExternalSorter::Record*> node;
    std::vector<ExternalSorter::Record*> rnode;// for creating root node 
    ExternalSorter::Record* r;

    while (true)
    {
	try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException) { break; }
	node.push_back(r);

	if (node.size() == b)
	{
	    Node* n = createNode(pTree, node, level);
	    node.clear();
	    pTree->writeNode(n);
	    rnode.push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, 0, 0));
	    // to standard output for detting those boundaries 
	    std::cout << n->m_identifier << " " << n->m_nodeMBR << std::endl;

	    delete n;
	}
    }

    if (! node.empty())
    {
	Node* n = createNode(pTree, node, level);
	pTree->writeNode(n);
	rnode.push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, 0, 0));
	// to standard output for detting those boundaries 
	std::cout << n->m_identifier << " " << n->m_nodeMBR << std::endl;
	delete n;
    }
    level++;

    // create final root node 
    Node* nr = createNode(pTree, rnode, level);
    pTree->writeNode(nr);
    pTree->m_rootID = nr->m_identifier;
    delete nr;

    pTree->m_stats.m_u32TreeHeight = level;
    pTree->storeHeader();
}


// BulkLoader
void BulkLoader::bulkLoadUsingRPLUS(
	SpatialIndex::RTree::RTree* pTree,
	IDataStream& stream,
	uint32_t partition_size,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
	) {
    if (! stream.hasNext())
	throw Tools::IllegalArgumentException(
		"RTree::BulkLoader::bulkLoadUsingRplus: Empty data stream given."
		);


    /* prepare the tree */
    NodePtr n = pTree->readNode(pTree->m_rootID);
    pTree->deleteNode(n.get());


//#ifndef NDEBUG
    std::cerr << "RTree::BulkLoader:R+ with K=" << partition_size << " , sorting data.."<< std::endl;
//#endif

    Tools::SmartPointer<ExternalSorter> ess = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
    Tools::SmartPointer<ExternalSorter> es = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(10000, 10000));

    const uint32_t DIM_X =0;
    const uint32_t DIM_Y =1;
    uint32_t dim = DIM_X; 

    while (stream.hasNext())
    {
	Data* d = reinterpret_cast<Data*>(stream.getNext());
	if (d == 0)
	    throw Tools::IllegalArgumentException(
		    "bulkLoadUsingRPLUS: RTree bulk load expects SpatialIndex::RTree::Data entries."
		    );

	es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, dim));
	d->m_pData = 0;
	delete d;
    }

    Region r = es->getUniverse();

// #ifndef NDEBUG
    std::cerr << "Spatial Universe: " << r << std::endl;
    std::cerr << "|collection| = " << es->getTotalEntries() << std::endl;
    std::cerr << "RTree::BulkLoader::R+ packing objects .." << std::endl;
// #endif

    std::vector<ExternalSorter::Record*> node;
    float cost [] = {0.0, 0.0};
    int iteration = 0; 
    id_type id = 0;
    
    while (true)
    {
	cost [0] = 0.0;
	cost [1] = 0.0;
	iteration++;

	if (es->getTotalEntries() <= partition_size) {

	    es->split(partition_size, dim, r,node);

// #ifndef NDEBUG
	    // last partition
	    std::cerr << "Iteration: " << iteration << "\tx-cost = " << cost [DIM_X] << "\ty-cost = " << cost [DIM_Y] << "\tRegion = " << r << std::endl;
	    std::cerr << "|collection| = " << es->getTotalEntries() << " , |partition| = " << node.size() << "." << std::endl;
// #endif
	    break; 
	}

	cost[DIM_X] = es->getCost(partition_size,DIM_X);
	cost[DIM_Y] = es->getCost(partition_size,DIM_Y);

	dim = (cost[DIM_X] <= cost[DIM_Y] )? DIM_X : DIM_Y ;
	es->split(partition_size, dim, r,node);

//#ifndef NDEBUG
	std::cerr << "Iteration: " << iteration << "\tx-cost = " << cost [DIM_X] << "\ty-cost = " << cost [DIM_Y] << "\tRegion = " << r << std::endl;
	std::cerr << "|collection| = " << es->getTotalEntries() << " , |partition| = " << node.size() << "." << std::endl;
//#endif

	ess->insert(new ExternalSorter::Record(r, id++, 0, 0, 0));
    } // end while 

    ess->insert(new ExternalSorter::Record(r, id++, 0, 0, 0));

    // STR style bulk loading for upper level nodes.
    pTree->m_stats.m_u64Data = ess->getTotalEntries();
    uint32_t level = 0;
    ess->sort();

    while (true)
    {
// #ifndef NDEBUG
	std::cerr << "RTree::BulkLoader::R+ Building level " << level << std::endl;
// #endif

	pTree->m_stats.m_nodesInLevel.push_back(0);

	Tools::SmartPointer<ExternalSorter> es2 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
	createLevel(pTree, ess, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
	ess = es2;

	if (ess->getTotalEntries() == 1) break;
	ess->sort();
    }

    pTree->m_stats.m_u32TreeHeight = level;
    pTree->storeHeader();
}


void BulkLoader::bulkLoadUsingSTR(
	SpatialIndex::RTree::RTree* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
	) {
    if (! stream.hasNext())
	throw Tools::IllegalArgumentException(
		"RTree::BulkLoader::bulkLoadUsingSTR: Empty data stream given."
		);

    NodePtr n = pTree->readNode(pTree->m_rootID);
    pTree->deleteNode(n.get());

#ifndef NDEBUG
    std::cerr << "RTree::BulkLoader: Sorting data." << std::endl;
#endif

    Tools::SmartPointer<ExternalSorter> es = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

    while (stream.hasNext())
    {
	Data* d = reinterpret_cast<Data*>(stream.getNext());
	if (d == 0)
	    throw Tools::IllegalArgumentException(
		    "bulkLoadUsingSTR: RTree bulk load expects SpatialIndex::RTree::Data entries."
		    );

	es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
	d->m_pData = 0;
	delete d;
    }
    es->sort();

    pTree->m_stats.m_u64Data = es->getTotalEntries();

    // create index levels.
    uint32_t level = 0;

    while (true)
    {
#ifndef NDEBUG
	std::cerr << "RTree::BulkLoader: Building level " << level << std::endl;
#endif

	pTree->m_stats.m_nodesInLevel.push_back(0);

	Tools::SmartPointer<ExternalSorter> es2 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
	createLevel(pTree, es, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
	es = es2;

	if (es->getTotalEntries() == 1) break;
	es->sort();
    }

    pTree->m_stats.m_u32TreeHeight = level;
    pTree->storeHeader();
}

void BulkLoader::createLevel(
	SpatialIndex::RTree::RTree* pTree,
	Tools::SmartPointer<ExternalSorter> es,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	Tools::SmartPointer<ExternalSorter> es2,
	uint32_t pageSize,
	uint32_t numberOfPages
	) {
    uint64_t b = (level == 0) ? bleaf : bindex;
    uint64_t P = static_cast<uint64_t>(std::ceil(static_cast<double>(es->getTotalEntries()) / static_cast<double>(b)));
    uint64_t S = static_cast<uint64_t>(std::ceil(std::sqrt(static_cast<double>(P))));

    if (S == 1 || dimension == pTree->m_dimension - 1 || S * b == es->getTotalEntries())
    {
	std::vector<ExternalSorter::Record*> node;
	ExternalSorter::Record* r;

	while (true)
	{
	    try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException) { break; }
	    node.push_back(r);

	    if (node.size() == b)
	    {
		Node* n = createNode(pTree, node, level);
		node.clear();
		pTree->writeNode(n);
		es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, 0, 0));
		pTree->m_rootID = n->m_identifier;
		// special case when the root has exactly bindex entries.
		delete n;
	    }
	}

	if (! node.empty())
	{
	    Node* n = createNode(pTree, node, level);
	    pTree->writeNode(n);
	    es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, 0, 0));
	    pTree->m_rootID = n->m_identifier;
	    delete n;
	}
    }
    else
    {
	bool bMore = true;

	while (bMore)
	{
	    ExternalSorter::Record* pR;
	    Tools::SmartPointer<ExternalSorter> es3 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

	    for (uint64_t i = 0; i < S * b; ++i)
	    {
		try { pR = es->getNextRecord(); }
		catch (Tools::EndOfStreamException) { bMore = false; break; }
		pR->m_s = dimension + 1;
		es3->insert(pR);
	    }
	    es3->sort();
	    createLevel(pTree, es3, dimension + 1, bleaf, bindex, level, es2, pageSize, numberOfPages);
	}
    }
}

Node* BulkLoader::createNode(SpatialIndex::RTree::RTree* pTree, std::vector<ExternalSorter::Record*>& e, uint32_t level)
{
    Node* n;

    if (level == 0) n = new Leaf(pTree, -1);
    else n = new Index(pTree, -1, level);

    for (size_t cChild = 0; cChild < e.size(); ++cChild)
    {
	n->insertEntry(e[cChild]->m_len, e[cChild]->m_pData, e[cChild]->m_r, e[cChild]->m_id);
	e[cChild]->m_pData = 0;
	delete e[cChild];
    }

    return n;
}
